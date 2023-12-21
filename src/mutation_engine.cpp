/*
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023 Alberto Casagrande <alberto.casagrande@uniud.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <filesystem>

#include <Rcpp.h>


#include <context_index.hpp>
#include <progress_bar.hpp>

#include "mutation_engine.hpp"

#include "genomic_data_storage.hpp"

struct MutationEngineSetup
{
  std::string description;
  std::filesystem::path directory;
  std::string reference_url;
  std::string SBS_url;
  size_t default_num_of_alleles;
  std::map<std::string, size_t> exceptions_on_allele_number;

  MutationEngineSetup(const std::string& description,
                  const std::filesystem::path& directory,
                  const std::string& reference_url,
                  const std::string& SBS_url,
                  const size_t& default_num_of_alleles,
                  const std::map<std::string, size_t>& exceptions_on_allele_number):
    description(description), directory(directory),
    reference_url(reference_url), SBS_url(SBS_url),
    default_num_of_alleles(default_num_of_alleles),
    exceptions_on_allele_number(exceptions_on_allele_number)
  {}
};

std::map<std::string, MutationEngineSetup> supported_setups{
  {
    "demo",
    {
      "A demostative set-up", "demo",
      "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2124/COSMIC_v3.4_SBS_GRCh38.txt",
      2, {{"X", 1}, {"Y", 1}}
    }
  },
  {
    "GRCh38",
    {
      "Homo sapiens (GRCh38)", "GRCh38",
      "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2124/COSMIC_v3.4_SBS_GRCh38.txt",
      2, {{"X", 1}, {"Y", 1}}
    }
  },
  {
    "GRCh37",
    {
      "Homo sapiens (GRCh37)", "GRCh37",
      "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh37.txt",
      2, {{"X", 1}, {"Y", 1}}
    }
  }
};

GenomicDataStorage setup_storage(const std::string& directory,
                                 const std::string& reference_url,
                                 const std::string& SBS_url)
{
  GenomicDataStorage storage(directory);

  storage.download_reference(reference_url);
  storage.download_SBS(SBS_url);

  storage.save_parameters();

  return storage;
}

GenomicDataStorage setup_storage(const std::string& setup_code)
{
  using namespace Rcpp;

  auto code_it = supported_setups.find(setup_code);
  if (code_it == supported_setups.end()) {
    std::ostringstream oss;
  
    oss << "\""+setup_code+"\" is an unknown code. "
        << "Supported codes are:" << std::endl;

    for (const auto& [s_code, setting]: supported_setups) {
      oss << " - \""+s_code+"\": " << setting.description << std::endl;
    }

    throw std::domain_error(oss.str());
  }

  return setup_storage(code_it->second.directory, 
                       code_it->second.reference_url,
                       code_it->second.SBS_url);
}

Rcpp::List MutationEngine::get_supported_setups()
{
  using namespace Rcpp;

  CharacterVector code_names(supported_setups.size()), code_descrs(supported_setups.size());

  size_t i=0;
  for (const auto& [s_code, setting]: supported_setups) {
    code_names[i] = s_code;
    code_descrs[i] = setting.description;
    ++i;
  }

  return DataFrame::create(_["name"]=code_names, _["description"]=code_descrs);
}

inline std::filesystem::path
get_context_index_path(const GenomicDataStorage& storage, const size_t context_sampling)
{
    return storage.get_directory()/std::string("context_index_"
                                                + std::to_string(context_sampling)
                                                + ".cif");
}

template<typename ABSOLUTE_GENOTYPE_POSITION = uint32_t>
Races::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION>
build_contex_index(const GenomicDataStorage& storage, const size_t context_sampling)
{
  using namespace Races;
  using namespace Races::Mutations;

  using Index = ContextIndex<ABSOLUTE_GENOTYPE_POSITION>;

  Index context_index;

  auto contex_index_filename = get_context_index_path(storage, context_sampling);

  if (std::filesystem::exists(contex_index_filename)) {
    Archive::Binary::In archive(contex_index_filename);
    archive.load(context_index, "context index");

    return context_index;
  }

  using namespace Rcpp;

  Rcout << "Building context index..." << std::endl << std::flush;

  std::vector<GenomicRegion> chr_regions;
  {
    UI::ProgressBar progress_bar;
    const auto reference_path = storage.get_reference_path();
    context_index = Index::build_index(reference_path, context_sampling,
                                       &progress_bar);
    chr_regions = context_index.get_chromosome_regions();
  }

  if (chr_regions.size() > 0) {
    Archive::Binary::Out archive(contex_index_filename);

    archive.save(context_index, "context index");
  }

  Rcout << "done" << std::endl;

  return context_index;
}

template<typename ABSOLUTE_GENOTYPE_POSITION>
std::map<Races::Mutations::ChromosomeId, size_t>
get_num_of_alleles(const Races::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION>& context_index,
                   const size_t& default_num_of_alleles,
                   const std::map<std::string, size_t>& alleles_num_exceptions)
{
  using namespace Races::Mutations;

  const auto chr_regions = context_index.get_chromosome_regions();

  std::map<ChromosomeId, size_t> alleles_per_chromosome;

  for (const auto& chr_region : chr_regions) {
    auto chr_id = chr_region.get_chromosome_id();

    alleles_per_chromosome[chr_id] = default_num_of_alleles;
  }

  for (const auto& [name, num_of_alleles] : alleles_num_exceptions) {
    auto chr_id = GenomicPosition::stochr(name);

    alleles_per_chromosome[chr_id] = num_of_alleles;
  }

  return alleles_per_chromosome;
}

std::map<std::string, Races::Mutations::MutationalSignature> 
load_SBS(const GenomicDataStorage& storage)
{
  std::ifstream is(storage.get_SBS_path());

  return Races::Mutations::MutationalSignature::read_from_stream(is);
}


void MutationEngine::init_mutation_engine(const GenomicDataStorage& storage,
                                          const size_t& default_num_of_alleles,
                                          const std::map<std::string, size_t>& alleles_num_exceptions,
                                          const size_t& context_sampling_rate)
{
  context_index = build_contex_index<MutationEngine::AbsGenotypePosition>(storage, context_sampling_rate);

  auto num_of_alleles = get_num_of_alleles(context_index, default_num_of_alleles, alleles_num_exceptions);

  auto SBS = load_SBS(storage);

  m_engine = Races::Mutations::MutationEngine(context_index, num_of_alleles, SBS);
}

MutationEngine::MutationEngine(const std::string& setup_name,
                               const size_t& context_sampling_rate)
{
  auto storage = setup_storage(setup_name);

  auto setup_cfg = supported_setups.at(setup_name);

  init_mutation_engine(storage, setup_cfg.default_num_of_alleles,
                       setup_cfg.exceptions_on_allele_number,
                       context_sampling_rate);
}

MutationEngine::MutationEngine(const std::string& directory,
                               const std::string& reference_url,
                               const std::string& SBS_url,
                               const size_t& default_num_of_alleles,
                               const std::map<std::string, size_t>& alleles_num_exceptions,
                               const size_t& context_sampling_rate)
{
  auto storage = setup_storage(directory, reference_url, SBS_url);

  init_mutation_engine(storage, default_num_of_alleles, alleles_num_exceptions, 
                       context_sampling_rate);
}

struct DummyTest
{
  template<typename RCPP_TYPE>
  void validate(const RCPP_TYPE& value) const
  {}
};

struct TestNonNegative
{
  template<typename RCPP_TYPE>
  void validate(const RCPP_TYPE& value)
  {
    auto c_value = Rcpp::as<double>(value);

    if (c_value<0) {
      throw std::runtime_error(std::to_string(c_value) + " should be non negative");
    }
  }
};

template<typename VALUE, typename TESTER=TestNonNegative>
std::map<std::string, VALUE>
get_map(const Rcpp::List& list)
{
  
  std::map<std::string, VALUE> c_map;
  if (list.size()==0) {
    return std::map<std::string, VALUE>();
  }

  TESTER tester;

  using namespace Rcpp;
  CharacterVector names = list.names();
  for (size_t i=0; i<list.size(); ++i) {
    tester.validate(list[i]);
    c_map[as<std::string>(names[i])] = as<VALUE>(list[i]);
  }

  return c_map;
}

MutationEngine 
MutationEngine::build_MutationEngine(const std::string& directory,
                                     const std::string& reference_url,
                                     const std::string& SBS_url,
                                     const size_t& default_num_of_alleles,
                                     const Rcpp::List& exceptions_on_allele_number,
                                     const std::string& setup_code,
                                     const size_t& context_sampling)
{
  if (setup_code!="") {
    if (directory!="" || reference_url!="" || SBS_url!="" || 
        default_num_of_alleles!=0 || exceptions_on_allele_number.size()!=0) {
      throw std::domain_error("when \"setup_code\" is provided, the parameters "
                              "\"directory\", \"reference_url\", \"SBS_url\", "
                              "\"default_num_of_alleles\", and "
                              "\"exceptions_on_allele_number\" must be empty");
    }

    return MutationEngine(setup_code, context_sampling);
  }

  if (directory=="" || reference_url=="" || SBS_url=="" || default_num_of_alleles==0) {
    throw std::domain_error("when \"setup_code\" is NOT provided, the parameters "
                            "\"directory\", \"reference_url\", \"SBS_url\", and "
                            "\"default_num_of_alleles\" are mandatory");
  }

  auto exceptions = get_map<size_t, TestNonNegative>(exceptions_on_allele_number);

  return MutationEngine(directory, reference_url, SBS_url, default_num_of_alleles,
                        exceptions, context_sampling);

}

void MutationEngine::add_coefficients(const double& time, const Rcpp::List& coefficients)
{
  auto m_coeffs = get_map<double, TestNonNegative>(coefficients);

  double sum{1};

  for (const auto& [sbs, coeff]: m_coeffs) {
    sum -= coeff;
  }

  if (std::abs(sum)>1e-13) {
    throw std::domain_error("The mutational coefficients must sum up to 1: "
                            "it sums up to "+std::to_string(sum));
  }

  m_engine.add(time, m_coeffs);
}

void MutationEngine::add_coefficients(const Rcpp::List& coefficients)
{
  add_coefficients(0, coefficients);
}


template<typename CPP_TYPE, typename RCPP_TYPE>
std::list<CPP_TYPE> get_ancestor_list(const Rcpp::List& rcpp_list)
{
  std::list<CPP_TYPE> cpp_list;

  for (size_t i=0; i<rcpp_list.size(); ++i) {
    cpp_list.push_back(static_cast<const CPP_TYPE&>(Rcpp::as<RCPP_TYPE>(rcpp_list[i])));
  }

  return cpp_list;
}

void MutationEngine::add_mutant(const std::string& mutant_name,
                                const Rcpp::List& species_rates,
                                const Rcpp::List& mutant_SNVs)
{
  Rcpp::List empty_list;

  add_mutant(mutant_name, species_rates, mutant_SNVs, empty_list);
}

void MutationEngine::add_mutant(const std::string& mutant_name,
                                const Rcpp::List& species_rates,
                                const Rcpp::List& mutant_SNVs,
                                const Rcpp::List& mutant_CNAs)
{
  auto c_snvs = get_ancestor_list<Races::Mutations::SNV, SNV>(mutant_SNVs);
  auto c_cnas = get_ancestor_list<Races::Mutations::CopyNumberAlteration, CNA>(mutant_CNAs);

  switch(species_rates.size()) {
    case 1:
      {
        double real_rate = species_rates[0];
        
        m_engine.add_mutant(mutant_name, {{"", real_rate}}, c_snvs, c_cnas);
      }
      break;
    case 2:
      {
        std::map<std::string, double> epi_rates;
        try {
          epi_rates = get_map<double, TestNonNegative>(species_rates);
        } catch (std::exception& e) {
          throw std::domain_error("The rates must be non-negative");
        };
        
        m_engine.add_mutant(mutant_name, epi_rates, c_snvs, c_cnas);
      }
      break;
    default:
      throw std::domain_error("Rates is either a real value or a named list "
                              "representing the mutation rates for the two "
                              "species.");
  }
}

template<typename OUT, typename ITERATOR>
OUT& show_list(OUT& os, ITERATOR it, ITERATOR last, const std::string& front="")
{
  for (; it != last; ++it) {
    os << front << *it << std::endl;
  }

  return os;
}

template<typename OUT, typename VALUE>
OUT& show_map(OUT& os, const std::map<std::string, VALUE>& c_map)
{
  os << "{";
  std::string sep = "";
  for (const auto& [key, value]: c_map) {
    os << sep << "\"" << key << "\": " << value;
    sep = ", ";
  }
  os << "}";

  return os;
}

void MutationEngine::show() const
{
  using namespace Rcpp;
  Rcout << "MutationEngine" << std::endl
        << " Species rates: ";
  
  const auto& m_properties = m_engine.get_mutational_properties();

  show_map(Rcout, m_properties.get_species_rates());

  Rcout << std::endl << std::endl << " Mutant genomic characterizations" << std::endl;
  for (const auto&[mutant_name, mutant_mutations]: m_properties.get_mutant_mutations()) {
    if (mutant_mutations.SNVs.size()>0 || mutant_mutations.CNAs.size()>0) {
      Rcout << "   \"" << mutant_name << "\"'s mutations " << std::endl;
    } else {
      Rcout << "   No mutations for \"" << mutant_name << "\"" << std::endl;
    }
    if (mutant_mutations.SNVs.size()>0) {
      Rcout << "     SNVs: " << std::endl;
            
      show_list(Rcout, mutant_mutations.SNVs.begin(), mutant_mutations.SNVs.end(), "       ");
    }
    if (mutant_mutations.CNAs.size()>0) {
      Rcout << "     CNAs: " << std::endl;

      show_list(Rcout, mutant_mutations.CNAs.begin(),  mutant_mutations.CNAs.end(), "       ");
    }
  }

  Rcout << std::endl << " Timed SBS coefficients" << std::endl;

  const auto& timed_mutational_coeffs = m_engine.get_timed_mutational_coefficients();
  auto coeffs_it = timed_mutational_coeffs.begin();
  while (coeffs_it != timed_mutational_coeffs.end()) {
    auto next_it = coeffs_it;
    ++next_it;
    if (next_it == timed_mutational_coeffs.end()) {
       Rcout << "   [" << coeffs_it->first << ", \u221E[: ";
    } else {
       Rcout << "   [" 
             << coeffs_it->first << ", " << next_it->first << "[: ";
    }
    show_map(Rcout, coeffs_it->second);
    Rcout << std::endl;

    coeffs_it = next_it;
  }
  Rcout << std::endl;
}
/*
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023-2024 Alberto Casagrande <alberto.casagrande@uniud.it>
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
#include <sstream>

#include <Rcpp.h>


#include <context_index.hpp>
#include <csv_reader.hpp>
#include <read_simulator.hpp>
#include <csv_reader.hpp>

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

Races::Mutations::GenomicRegion get_CNA_region(const Races::IO::CSVReader::CSVRow& row, const size_t& row_num)
{
  using namespace Races::Mutations;

  ChromosomeId chr_id;    
  try {
    chr_id = GenomicPosition::stochr(row.get_field(0).substr(3));
  } catch (std::invalid_argument) {
    throw std::domain_error("Unknown chromosome specification " + row.get_field(1) 
                            + " in row number " + std::to_string(row_num) 
                            + ".");
  }

  uint32_t begin_pos;         
  try {
    begin_pos = stoul(row.get_field(1));
  } catch (std::invalid_argument) {
    throw std::domain_error("Unknown begin specification " + row.get_field(1) 
                            + " in row number " + std::to_string(row_num) 
                            + ".");
  }

  GenomicPosition pos(chr_id, begin_pos);

  uint32_t end_pos;                
  try {
    end_pos = stoul(row.get_field(2));
  } catch (std::invalid_argument) {
    throw std::domain_error("Unknown end specification " + row.get_field(2) 
                            + " in row number " + std::to_string(row_num) 
                            + ".");
  }

  if (begin_pos>end_pos) {
    throw std::domain_error("The CNA begin lays after the end in row number "
                            + std::to_string(row_num));
  }
  
  return {pos, end_pos+1-begin_pos};
}


std::vector<Races::Mutations::CopyNumberAlteration> load_passenger_CNAs(const std::filesystem::path& CNAs_csv,
                                                                        const std::string& tumor_type)
{
  using namespace Races::Mutations;

  std::vector<CopyNumberAlteration> CNAs;

  Races::IO::CSVReader csv_reader(CNAs_csv);

  size_t row_num{2};
  for (const auto& row : csv_reader) {
    if (row.size()<6) {
      throw std::runtime_error("The CNA CSV must contains at least 6 columns");
    }
    if ((tumor_type=="") || (row.get_field(5) == tumor_type)) {
      const auto region = get_CNA_region(row, row_num);

      const auto major = row.get_field(3);
      try {
        if (major=="NA" || (stoi(major)>1)) {
          CNAs.push_back({region, CopyNumberAlteration::Type::AMPLIFICATION});
        }
      } catch (std::invalid_argument) {
        throw std::domain_error("Unknown major specification " + major 
                                + " in row number " + std::to_string(row_num) 
                                + ".");
      }

      const auto minor = row.get_field(4);
      try {
        if (minor=="NA" || (stoi(minor)<1)) {
          CNAs.push_back({region, CopyNumberAlteration::Type::DELETION});
        }
      } catch (std::invalid_argument) {
        throw std::domain_error("Unknown minor specification " + major 
                                + " in row number " + std::to_string(row_num) 
                                + ".");
      }
    }

    ++row_num;
  }

  return CNAs;
}

std::filesystem::path get_passenger_data_file_path()
{
  using namespace Rcpp;

  Function system_file("system.file");

  auto CNA_file = system_file("extdata", "passenger_CNAs.csv", _["package"]="rRACES");

  return as<std::string>(CNA_file);
}

std::list<std::string> MutationEngine::get_supported_tumor_types()
{
  auto CNA_file = get_passenger_data_file_path();

  std::set<std::string> tumor_types;

  Races::IO::CSVReader csv_reader(CNA_file);
  for (const auto& row : csv_reader) {
    tumor_types.insert(row.get_field(5));
  }

  return {tumor_types.begin(),tumor_types.end()};
}

std::vector<Races::Mutations::CopyNumberAlteration>
load_passenger_CNAs(const std::string& tumor_type="")
{
  auto CNA_file = get_passenger_data_file_path();

  return load_passenger_CNAs(CNA_file, tumor_type);
}

void MutationEngine::init_mutation_engine(const GenomicDataStorage& storage,
                                          const size_t& default_num_of_alleles,
                                          const std::map<std::string, size_t>& alleles_num_exceptions,
                                          const size_t& context_sampling_rate,
                                          const std::string& tumor_type)
{
  reference_path = storage.get_reference_path();
  SBS_path = storage.get_SBS_path();

  context_index = build_contex_index<MutationEngine::AbsGenotypePosition>(storage, context_sampling_rate);

  auto num_of_alleles = get_num_of_alleles(context_index, default_num_of_alleles, alleles_num_exceptions);

  auto SBS = load_SBS(storage);

  auto passenger_CNAs = load_passenger_CNAs(tumor_type);

  m_engine = Races::Mutations::MutationEngine(context_index, num_of_alleles, SBS, passenger_CNAs);
}

MutationEngine::MutationEngine(const std::string& setup_name,
                               const size_t& context_sampling_rate,
                               const std::string& tumor_type)
{
  auto storage = setup_storage(setup_name);

  auto setup_cfg = supported_setups.at(setup_name);

  init_mutation_engine(storage, setup_cfg.default_num_of_alleles,
                       setup_cfg.exceptions_on_allele_number,
                       context_sampling_rate, tumor_type);
}

MutationEngine::MutationEngine(const std::string& directory,
                               const std::string& reference_url,
                               const std::string& SBS_url,
                               const size_t& default_num_of_alleles,
                               const std::map<std::string, size_t>& alleles_num_exceptions,
                               const size_t& context_sampling_rate,
                               const std::string& tumor_type)
{
  auto storage = setup_storage(directory, reference_url, SBS_url);

  init_mutation_engine(storage, default_num_of_alleles, alleles_num_exceptions, 
                       context_sampling_rate, tumor_type);
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
                                     const size_t& context_sampling,
                                     const std::string& tumor_type)
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
                        exceptions, context_sampling, tumor_type);

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

void MutationEngine::add_exposure(const double& time, const Rcpp::List& exposure)
{
  auto c_exposure = get_map<double, TestNonNegative>(exposure);

  double sum{1};

  for (const auto& [sbs, coeff]: c_exposure) {
    sum -= coeff;
  }

  if (std::abs(sum)>1e-13) {
    std::ostringstream oss;

    oss << "The exposure must sum up to 1: ";
    show_map(oss, c_exposure);
    oss << " sums up to " << sum;

    throw std::domain_error(oss.str());
  }

  m_engine.add(time, c_exposure);
}

void MutationEngine::add_exposure(const Rcpp::List& exposure)
{
  add_exposure(0, exposure);
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
                                const Rcpp::List& epistate_passenger_rates,
                                const Rcpp::List& driver_SNVs)
{
  Rcpp::List empty_list;

  add_mutant(mutant_name, epistate_passenger_rates, driver_SNVs, empty_list);
}

double get_non_negative(const Rcpp::List& values,
                        const std::string& message=": expected non-negative value")
{
  if (values.size()!=1) {
    throw std::runtime_error("Expected one non-negative value. Got a list of "
                             + std::to_string(values.size()) + " values.");
  }
  auto c_value = Rcpp::as<double>(values[0]);

  if (c_value<0) {
    throw std::runtime_error(std::to_string(c_value) + message);
  }

  return c_value;
}

bool contains_passenger_rates(const Rcpp::List& list)
{
  Rcpp::CharacterVector names = list.names();
  for (size_t i=0; i<list.size(); ++i) {
    if (names[i]!="SNV" && names[i]!="CNA") {
      return false;
    }
  }

  return true;
}

Races::Mutations::PassengerRates
get_passenger_rates(const Rcpp::List& list)
{
  Races::Mutations::PassengerRates p_rates;

  if (!list.hasAttribute("names")) {
    throw std::runtime_error("Passenger rates list must be a named list whose names "
                             "are in the set {\"SNV\", \"CNA\"}.");
  }

  Rcpp::CharacterVector names = list.names();
  for (size_t i=0; i<list.size(); ++i) {
    if (names[i]=="SNV") {
      p_rates.snv = get_non_negative(list[i], ": SNV rates must be non-negative");
    } else {
      if (names[i]=="CNA") {
        p_rates.cna = get_non_negative(list[i], ": CNA rates must be non-negative");
      } else {
        throw std::runtime_error("\"" + names[i] + "\" is an unsupported name in "
                                 + "passenger rates list.");
      }
    }
  }

  return p_rates;
}

std::map<std::string, Races::Mutations::PassengerRates>
get_epistate_passenger_rates(const Rcpp::List& list)
{
  std::map<std::string, Races::Mutations::PassengerRates> ep_rates;

  if (!list.hasAttribute("names")) {
    throw std::runtime_error("Epistate passenger rates list must be a "
                             "named list whose names are epistates, "
                             "i.e., either \"+\" or \"-\".");
  }

  Rcpp::CharacterVector names = list.names();
  for (size_t i=0; i<list.size(); ++i) {
    if (names[i]=="+") {
      ep_rates["+"] = get_passenger_rates(list[i]);
    } else {
      if (names[i]=="-") {
        ep_rates["-"] = get_passenger_rates(list[i]);
      } else {
        throw std::runtime_error("\"" + names[i] + "\" is an unsupported name in "
                                 + "epistate passenger rate list.");
      }
    }
  }

  return ep_rates;
}

struct FilterNonChromosomeSequence : public Races::IO::FASTA::SequenceFilter
{
    Races::Mutations::ChromosomeId last_chr_id;

    inline bool operator()(const std::string& header)
    {
        return !Races::IO::FASTA::is_chromosome_header(header, last_chr_id);
    } 
};

using SNV_iterator = std::list<std::list<Races::Mutations::SNV>::iterator>;

void check_wrong_chromosome_SNV(const std::map<Races::Mutations::ChromosomeId, SNV_iterator>& SNV_partition)
{
  if (SNV_partition.size()>0) {
    std::ostringstream oss;

    std::string sep="";
    size_t counter{0};
    for (const auto& [chr, SNV_class] : SNV_partition) {
      for (auto& SNV_it : SNV_class) {
        oss << sep << *SNV_it;

        if (sep.size()==0) {
          sep = ", ";
        }
        ++counter;
      }
    }

    throw std::runtime_error((counter>1?"SNVs ":"SNV ")
                            + oss.str()
                            + (counter>1?" belong":" belongs")
                            + (SNV_partition.size()>1?" to unknown chromosomes":
                                                      " to an unknown chromosome"));
  }
}

void check_SNVs_and_retrieve_missing_contexts(const std::filesystem::path& fasta_filename, std::list<Races::Mutations::SNV>& SNVs)
{
  std::map<Races::Mutations::ChromosomeId, SNV_iterator> SNV_partition;

  for (auto it=SNVs.begin(); it != SNVs.end(); ++it) {
    SNV_partition[it->chr_id].push_back(it);
  }

  Races::IO::FASTA::Sequence chr_seq;

  FilterNonChromosomeSequence filter;

  size_t SNV_to_check{SNVs.size()};
  std::ifstream fasta_stream(fasta_filename);
  while (Races::IO::FASTA::Sequence::read(fasta_stream, chr_seq, filter)) {
    auto chr_id = filter.last_chr_id;

    auto found = SNV_partition.find(chr_id);
    if (found != SNV_partition.end()) {
      for (auto& SNV_it : found->second) {
        if (SNV_it->position >= chr_seq.nucleotides.size()) {
          std::ostringstream oss;

          oss << "The SNV context of " << *SNV_it
              << " does not lay into the chromosome." << std::endl;
          throw std::out_of_range(oss.str());
        }

        auto candidate_context = chr_seq.nucleotides.substr(SNV_it->position-2, 3);
        if (!SNV_it->context.is_defined()) {
          if (candidate_context.find('N') == std::string::npos) {
            SNV_it->context = candidate_context;
          }
        } else {
          if (SNV_it->context.get_sequence() != candidate_context) {
            std::ostringstream oss;

            oss << "Reference contains sequence \""
                << candidate_context << "\" instead of \"" 
                << SNV_it->context.get_sequence()
                << "\" as context for " 
                << *SNV_it << "." << std::endl;
            throw std::out_of_range(oss.str());
          }
        }
        --SNV_to_check;
      }
      if (SNV_to_check == 0) {
        return;
      }

      SNV_partition.erase(found);
    }
  }

  check_wrong_chromosome_SNV(SNV_partition);
}

void MutationEngine::add_mutant(const std::string& mutant_name,
                                const Rcpp::List& epistate_passenger_rates,
                                const Rcpp::List& driver_SNVs,
                                const Rcpp::List& driver_CNAs)
{
  auto c_snvs = get_ancestor_list<Races::Mutations::SNV, SNV>(driver_SNVs);
  auto c_cnas = get_ancestor_list<Races::Mutations::CopyNumberAlteration, CNA>(driver_CNAs);

  check_SNVs_and_retrieve_missing_contexts(reference_path, c_snvs);

  if (contains_passenger_rates(epistate_passenger_rates)) {
    auto p_rates = get_passenger_rates(epistate_passenger_rates);
    m_engine.add_mutant(mutant_name, {{"", p_rates}}, c_snvs, c_cnas);

    return;
  }

  auto epi_rates = get_epistate_passenger_rates(epistate_passenger_rates);
  m_engine.add_mutant(mutant_name, epi_rates, c_snvs, c_cnas);
}

PhylogeneticForest MutationEngine::place_mutations(const SamplesForest& forest, 
                                                   const size_t& num_of_preneoplatic_mutations,
                                                   const int seed)
{
  Races::UI::ProgressBar progress_bar;

  progress_bar.set_message("Placing mutations");

  auto phylo_forest = m_engine.place_mutations(forest, num_of_preneoplatic_mutations, progress_bar, seed);

  progress_bar.set_message("Mutations placed");

  return {std::move(phylo_forest), reference_path, m_engine.get_timed_exposures()};
}

Rcpp::List MutationEngine::get_SBS_dataframe()
{
  Rcpp::Function read_delim("read.delim");
  
  std::cout << SBS_path << std::endl;

  return read_delim(std::string(SBS_path), Rcpp::_["quote"]="");
}

template<typename OUT, typename ITERATOR>
OUT& show_list(OUT& os, ITERATOR it, ITERATOR last, const std::string& front="")
{
  for (; it != last; ++it) {
    os << front << *it << std::endl;
  }

  return os;
}

void MutationEngine::show() const
{
  using namespace Rcpp;
  Rcout << "MutationEngine" << std::endl
        << " Passenger rates";
  
  const auto& m_properties = m_engine.get_mutational_properties();

  for (const auto& [species_name, p_rates] : m_properties.get_passenger_rates()) {
    Rcout << std::endl << "   \"" << species_name << "\": {";
    std::string sep;
    if (p_rates.snv>0) {
      Rcout << "SNV: " << p_rates.snv;
      sep = ", ";
    }
    if (p_rates.cna>0) {
      Rcout << sep << "CNA: " << p_rates.cna;
    }
    Rcout << "}";
  }

  Rcout << std::endl << std::endl << " Driver mutations" << std::endl;
  for (const auto&[mutant_name, driver_mutations]: m_properties.get_driver_mutations()) {
    if (driver_mutations.SNVs.size()>0 || driver_mutations.CNAs.size()>0) {
      if (driver_mutations.SNVs.size()>0) {
        Rcout << "   \"" << mutant_name << "\" SNVs: " << std::endl;
              
        show_list(Rcout, driver_mutations.SNVs.begin(), driver_mutations.SNVs.end(), "       ");
      }
      if (driver_mutations.CNAs.size()>0) {
        Rcout << "   \"" << mutant_name << "\" CNAs: " << std::endl;

        show_list(Rcout, driver_mutations.CNAs.begin(),  driver_mutations.CNAs.end(), "       ");
      }
    } else {
      Rcout << "   No driver mutations for \"" << mutant_name << "\"" << std::endl;
    }
  }

  Rcout << std::endl << " Timed Exposure" << std::endl;

  const auto& timed_exposures = m_engine.get_timed_exposures();
  auto coeffs_it = timed_exposures.begin();
  while (coeffs_it != timed_exposures.end()) {
    auto next_it = coeffs_it;
    ++next_it;
    if (next_it == timed_exposures.end()) {
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
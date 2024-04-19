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

#include <utils.hpp>

#include "mutation_engine.hpp"

#include "genomic_data_storage.hpp"

struct MutationEngineSetup
{
  std::string description;
  std::filesystem::path directory;
  std::string reference_url;
  std::string SBS_url;
  std::string drivers_url;
  std::string passenger_CNAs_url;
  std::string germline_url;

  MutationEngineSetup(const std::string& description,
                  const std::filesystem::path& directory,
                  const std::string& reference_url,
                  const std::string& SBS_url,
                  const std::string& drivers_url,
                  const std::string& passenger_CNAs_url,
                  const std::string& germline_url):
    description(description), directory(directory), reference_url(reference_url), 
    SBS_url(SBS_url), drivers_url(drivers_url), 
    passenger_CNAs_url(passenger_CNAs_url), germline_url(germline_url)
  {}
};

std::map<std::string, MutationEngineSetup> supported_setups{
  {
    "demo",
    {
      "A demostative set-up", "demo",
      "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh37.txt",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/driver_mutations_hg19.csv",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/passenger_CNAs_hg19.csv",
      "https://www.dropbox.com/scl/fi/g9oloxkip18tr1rm6wjve/germline_data_demo.tar.gz?rlkey=15jshuld3bqgyfcs7fa0bzqeo&dl=1"
    }
  },
  {
    "GRCh38",
    {
      "Homo sapiens (GRCh38)", "GRCh38",
      "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2124/COSMIC_v3.4_SBS_GRCh38.txt",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/driver_mutations_hg38.csv",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/passenger_CNAs_hg38.csv",
      "https://www.dropbox.com/scl/fi/3rs2put4wde3objxmhvjc/germline_data_hg38.tar.gz?rlkey=imawitklf8d6zphz9ugriv4qm&dl=1"
    }
  },
  {
    "GRCh37",
    {
      "Homo sapiens (GRCh37)", "GRCh37",
      "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh37.txt",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/driver_mutations_hg19.csv",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/passenger_CNAs_hg19.csv",
      "https://www.dropbox.com/scl/fi/ckj7k0db0v1qf0o8la2yx/germline_data_hg19.tar.gz?rlkey=aanaz7n9v1bvmfvuqvqamc76o&dl=1"
    }
  }
};

GenomicDataStorage setup_storage(const std::string& directory,
                                 const std::string& reference_source,
                                 const std::string& SBS_source,
                                 const std::string& drivers_source,
                                 const std::string& passengers_CNA_source,
                                 const std::string& germline_source)
{
  GenomicDataStorage storage(directory, reference_source, SBS_source,
                             drivers_source, passengers_CNA_source, 
                             germline_source);

  storage.save_sources();

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

  auto directory = to_string(code_it->second.directory);

  return setup_storage(directory,
                       code_it->second.reference_url,
                       code_it->second.SBS_url,
                       code_it->second.drivers_url,
                       code_it->second.passenger_CNAs_url,
                       code_it->second.germline_url);
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
    archive.load(context_index, "context index", Rcpp::Rcout);

    return context_index;
  }

  using namespace Rcpp;

  Rcout << "Building context index..." << std::endl << std::flush;

  std::set<GenomicRegion> regions_to_avoid;

  auto drivers_path = storage.get_driver_mutations_path();
  if (std::filesystem::exists(drivers_path)) {
    auto driver_storage = DriverStorage::load(drivers_path);

    for (const auto& [name, snv] : driver_storage.get_SNVs()) {
        regions_to_avoid.emplace(snv, 1);
    }
  }

  std::list<GenomicRegion> chr_regions;
  {
    UI::ProgressBar progress_bar(Rcpp::Rcout);
    const auto reference_path = storage.get_reference_path();
    context_index = Index::build_index(reference_path, regions_to_avoid,
                                       context_sampling, &progress_bar);
    chr_regions = context_index.get_chromosome_regions();
  }

  if (chr_regions.size() > 0) {
    Archive::Binary::Out archive(contex_index_filename);

    archive.save(context_index, "context index", Rcpp::Rcout);
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
  } catch (std::invalid_argument const&) {
    throw std::domain_error("Unknown chromosome specification " + row.get_field(1) 
                            + " in row number " + std::to_string(row_num) 
                            + ".");
  }

  uint32_t begin_pos;         
  try {
    begin_pos = stoul(row.get_field(1));
  } catch (std::invalid_argument const&) {
    throw std::domain_error("Unknown begin specification " + row.get_field(1) 
                            + " in row number " + std::to_string(row_num) 
                            + ".");
  }

  GenomicPosition pos(chr_id, begin_pos);

  uint32_t end_pos;                
  try {
    end_pos = stoul(row.get_field(2));
  } catch (std::invalid_argument const&) {
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


std::vector<Races::Mutations::CNA> load_passenger_CNAs(const std::filesystem::path& CNAs_csv,
                                                                        const std::string& tumor_type)
{
  std::vector<Races::Mutations::CNA> CNAs;

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
          CNAs.emplace_back(region.get_begin(), region.size(),
                            CNA::Type::AMPLIFICATION);
        }
      } catch (std::invalid_argument const&) {
        throw std::domain_error("Unknown major specification " + major 
                                + " in row number " + std::to_string(row_num) 
                                + ".");
      }

      const auto minor = row.get_field(4);
      try {
        if (minor=="NA" || (stoi(minor)<1)) {
          CNAs.emplace_back(region.get_begin(), region.size(),
                            CNA::Type::DELETION);
        }
      } catch (std::invalid_argument const&) {
        throw std::domain_error("Unknown minor specification " + major 
                                + " in row number " + std::to_string(row_num) 
                                + ".");
      }
    }

    ++row_num;
  }

  return CNAs;
}

void MutationEngine::init_mutation_engine()
{
  context_index = build_contex_index<MutationEngine::AbsGenotypePosition>(storage, context_sampling);

  reset();
}

MutationEngine::MutationEngine(const std::string& setup_name,
                               const std::string& germline_subject,
                               const size_t& context_sampling,
                               const std::string& tumor_type):
  storage(setup_storage(setup_name)), germline_subject(germline_subject), 
  context_sampling(context_sampling), tumor_type(tumor_type)
{
  auto setup_cfg = supported_setups.at(setup_name);

  init_mutation_engine();
}

MutationEngine::MutationEngine(const std::string& directory,
                               const std::string& reference_source,
                               const std::string& SBS_source,
                               const std::string& drivers_source,
                               const std::string& passenger_CNAs_source,
                               const std::string& germline_source,
                               const std::string& germline_subject,
                               const size_t& context_sampling,
                               const std::string& tumor_type):
  storage(setup_storage(directory, reference_source, SBS_source, 
                        drivers_source, passenger_CNAs_source, 
                        germline_source)),
  germline_subject(germline_subject), context_sampling(context_sampling),
  tumor_type(tumor_type)
{
  init_mutation_engine();
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
  const size_t list_size = static_cast<size_t>(list.size());
  for (size_t i=0; i<list_size; ++i) {
    tester.validate(list[i]);
    c_map[as<std::string>(names[i])] = as<VALUE>(list[i]);
  }

  return c_map;
}

MutationEngine 
MutationEngine::build_MutationEngine(const std::string& directory,
                                     const std::string& reference_source,
                                     const std::string& SBS_source,
                                     const std::string& drivers_source,
                                     const std::string& passenger_CNAs_source,
                                     const std::string& germline_source,
                                     const std::string& setup_code,
                                     const std::string& germline_subject,
                                     const size_t& context_sampling,
                                     const std::string& tumor_type)
{
  if (setup_code!="") {
    if (directory!="" || reference_source!="" || SBS_source!="" 
         || drivers_source!="" || passenger_CNAs_source !="" 
         || germline_source !="") {
      throw std::domain_error("when \"setup_code\" is provided, the parameters "
                              "\"directory\", \"reference_src\", \"SBS_src\", "
                              "\"drivers_src\", \"passenger_CNAs_src\" and "
                              "\"germline_src\" must be avoided.");
    }

    return MutationEngine(setup_code, germline_subject,
                          context_sampling, tumor_type);
  }

  if (directory=="" || reference_source=="" || SBS_source=="" 
      || passenger_CNAs_source== "" || germline_source== "") {
    throw std::domain_error("when \"setup_code\" is NOT provided, the parameters "
                            "\"directory\", \"reference_src\", \"SBS_src\", "
                            "\"passenger_CNAs_src\", and \"germline_src\" are "
                            "mandatory.");
  }

  return MutationEngine(directory, reference_source, SBS_source, drivers_source,
                        passenger_CNAs_source, germline_source, germline_subject, 
                        context_sampling, tumor_type);

}

template<typename VALUE>
std::ostream& show_map(std::ostream& os, const std::map<std::string, VALUE>& c_map)
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
std::list<CPP_TYPE> get_super_object_list(const Rcpp::List& rcpp_list)
{
  std::list<CPP_TYPE> cpp_list;

  const size_t list_size = static_cast<size_t>(rcpp_list.size());
  for (size_t i=0; i<list_size; ++i) {
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
  const size_t list_size = static_cast<size_t>(list.size());
  for (size_t i=0; i<list_size; ++i) {
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
  const size_t list_size = static_cast<size_t>(list.size());
  for (size_t i=0; i<list_size; ++i) {
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
  const size_t list_size = static_cast<size_t>(list.size());
  for (size_t i=0; i<list_size; ++i) {
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

using SNVSpec = Races::Mutations::MutationSpec<Races::Mutations::SNV>;
using SNV_iterator = std::list<std::list<SNVSpec>::iterator>;

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

inline std::ifstream::pos_type filesize(const std::filesystem::path& fasta_filename)
{
    std::ifstream in(fasta_filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

void retrieve_missing_references(const std::string& mutant_name, 
                                 const std::filesystem::path& fasta_filename,
                                 std::list<SNVSpec>& SNVs)
{
  Races::UI::ProgressBar progress_bar(Rcpp::Rcout);

  std::map<Races::Mutations::ChromosomeId, SNV_iterator> SNV_partition;

  size_t SNV_to_check{0};
  for (auto it=SNVs.begin(); it != SNVs.end(); ++it) {
    if (it->ref_base == '?') {
      SNV_partition[it->chr_id].push_back(it);
      ++SNV_to_check;
    }
  }

  Races::IO::FASTA::Sequence chr_seq;

  FilterNonChromosomeSequence filter;

  progress_bar.set_message("Retrieving \""+ mutant_name + "\" SNVs");

  const auto fasta_size = filesize(fasta_filename);
  std::ifstream fasta_stream(fasta_filename);
  while (SNV_to_check>0 && Races::IO::FASTA::Sequence::read(fasta_stream, chr_seq, filter, progress_bar)) {
    auto chr_id = filter.last_chr_id;

    progress_bar.set_progress((100*fasta_stream.tellg())/fasta_size);

    auto found = SNV_partition.find(chr_id);
    if (found != SNV_partition.end()) {
      for (auto& SNV_it : found->second) {
        if (SNV_it->position >= chr_seq.nucleotides.size()) {
          std::ostringstream oss;

          oss << "The SNV context of " << *SNV_it
              << " does not lay into the chromosome." << std::endl;
          throw std::out_of_range(oss.str());
        }

        const auto& candidate_ref = chr_seq.nucleotides[SNV_it->position-1];
        if (SNV_it->ref_base == '?') {
          SNV_it->ref_base = candidate_ref;
        }
        --SNV_to_check;
      }

      SNV_partition.erase(found);
    }
  }

  check_wrong_chromosome_SNV(SNV_partition);

  progress_bar.set_progress(100, "\"" + mutant_name + "\" SNVs retrieved");
}

void MutationEngine::add_mutant(const std::string& mutant_name,
                                const Rcpp::List& epistate_passenger_rates,
                                const Rcpp::List& driver_SNVs,
                                const Rcpp::List& driver_CNAs)
{
  auto c_snvs = get_super_object_list<SNVSpec, SNV>(driver_SNVs);
  auto c_cnas = get_super_object_list<Races::Mutations::CNA, CNA>(driver_CNAs);

  retrieve_missing_references(mutant_name, storage.get_reference_path(), c_snvs);

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
  Races::UI::ProgressBar progress_bar(Rcpp::Rcout);

  progress_bar.set_message("Placing mutations");

  auto phylo_forest = m_engine.place_mutations(forest, num_of_preneoplatic_mutations, progress_bar, seed);

  progress_bar.set_message("Mutations placed");

  return {std::move(phylo_forest), germline_subject, storage.get_reference_path(),
          m_engine.get_timed_exposures()};
}

Rcpp::List MutationEngine::get_SBS_dataframe()
{
  Rcpp::Function read_delim("read.delim");

  return read_delim(to_string(storage.get_SBS_path()), Rcpp::_["quote"]="");
}

template<typename ITERATOR>
std::ostream& show_list(std::ostream& os, ITERATOR it, ITERATOR last, const std::string& front="")
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

void MutationEngine::rebuild_context_index()
{
  auto context_index_path = get_context_index_path(storage, context_sampling);

  std::filesystem::remove(context_index_path);

  context_index = build_contex_index(storage, context_sampling);
}

void MutationEngine::set_context_sampling(const size_t& context_sampling)
{
  this->context_sampling = context_sampling;

  context_index = build_contex_index(storage, context_sampling);
}

void MutationEngine::reset(const bool full)
{
  using namespace Races::Mutations;

  MutationalProperties mutational_properties;
  std::map<Races::Time, Exposure> timed_exposures;

  if (!full) {
    mutational_properties = m_engine.get_mutational_properties();
    timed_exposures = m_engine.get_timed_exposures();
  }

  auto SBS = load_SBS(storage);

  auto passenger_CNAs = load_passenger_CNAs(storage.get_passenger_CNAs_path(),
                                            tumor_type);

  auto driver_storage = DriverStorage::load(storage.get_driver_mutations_path());

  auto& germline_storage = storage.get_germline_storage();

  if (germline_subject == "") {
    auto germline_subjects = germline_storage.get_population();

    if (germline_subjects.size() == 0) {
      throw std::runtime_error("No germline subject available.");
    }

    germline_subject = germline_subjects[0].name;
  }

  auto germline = germline_storage.get_germline(germline_subject);

  m_engine = Races::Mutations::MutationEngine(context_index, SBS,
                                              mutational_properties, germline,
                                              driver_storage, passenger_CNAs);

  for (const auto& [time, exposure] : timed_exposures) {
    m_engine.add(time, exposure);
  }
}

void MutationEngine::set_germline_subject(const std::string& germline_subject)
{
  storage.get_germline_storage().get_subject(germline_subject);

  this->germline_subject = germline_subject;

  reset(false);
}

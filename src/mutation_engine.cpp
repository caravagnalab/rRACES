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

#include "utility.hpp"

#include "genomic_data_storage.hpp"

using SIDSpec = RACES::Mutations::MutationSpec<RACES::Mutations::SID>;
using SID_iterator = std::list<std::list<SIDSpec>::iterator>;

struct MutationEngineSetup
{
  std::string description;
  std::filesystem::path directory;
  std::string reference_url;
  std::string SBS_signatures_url;
  std::string indel_signatures_url;
  std::string drivers_url;
  std::string passenger_CNAs_url;
  std::string germline_url;

  MutationEngineSetup(const std::string& description,
                      const std::filesystem::path& directory,
                      const std::string& reference_url,
                      const std::string& SBS_signatures_url,
                      const std::string& indel_signatures_url,
                      const std::string& drivers_url,
                      const std::string& passenger_CNAs_url,
                      const std::string& germline_url):
    description(description), directory(directory), reference_url(reference_url),
    SBS_signatures_url(SBS_signatures_url),
    indel_signatures_url(indel_signatures_url),
    drivers_url(drivers_url),
    passenger_CNAs_url(passenger_CNAs_url), germline_url(germline_url)
  {}
};

std::map<std::string, MutationEngineSetup> supported_setups{
  {
    "demo",
    {
      "A demostative set-up", "demo",
      "https://ftp.ensembl.org/pub/grch37/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh37.txt",
      "https://cancer.sanger.ac.uk/signatures/documents/2121/COSMIC_v3.4_ID_GRCh37.txt",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/driver_mutations_hg19.csv",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/passenger_CNAs_hg19.csv",
      "https://zenodo.org/records/13166780/files/germline_data_demo.tar.gz"
    }
  },
  {
    "GRCh38",
    {
      "Homo sapiens (GRCh38)", "GRCh38",
      "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2124/COSMIC_v3.4_SBS_GRCh38.txt",
      "https://cancer.sanger.ac.uk/signatures/documents/2121/COSMIC_v3.4_ID_GRCh37.txt",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/driver_mutations_hg38.csv",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/passenger_CNAs_hg38.csv",
      "https://zenodo.org/records/13166780/files/germline_data_hg38.tar.gz"
    }
  },
  {
    "GRCh37",
    {
      "Homo sapiens (GRCh37)", "GRCh37",
      "https://ftp.ensembl.org/pub/grch37/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh37.txt",
      "https://cancer.sanger.ac.uk/signatures/documents/2121/COSMIC_v3.4_ID_GRCh37.txt",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/driver_mutations_hg19.csv",
      "https://raw.githubusercontent.com/caravagnalab/rRACES/main/inst/extdata/passenger_CNAs_hg19.csv",
      "https://zenodo.org/records/13166780/files/germline_data_hg19.tar.gz"
    }
  }
};

GenomicDataStorage setup_storage(const std::string& directory,
                                 const std::string& reference_source,
                                 const std::string& SBS_signatures_source,
                                 const std::string& indel_signatures_source,
                                 const std::string& drivers_source,
                                 const std::string& passengers_CNA_source,
                                 const std::string& germline_source)
{
  GenomicDataStorage storage(directory, reference_source, SBS_signatures_source,
                             indel_signatures_source,
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
                       code_it->second.SBS_signatures_url,
                       code_it->second.indel_signatures_url,
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
RACES::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION>
build_contex_index(const GenomicDataStorage& storage, const size_t context_sampling,
                   const bool& quiet)
{
  using namespace RACES;
  using namespace RACES::Mutations;

  using Index = ContextIndex<ABSOLUTE_GENOTYPE_POSITION>;

  Index context_index;

  auto contex_index_filename = get_context_index_path(storage, context_sampling);

  if (std::filesystem::exists(contex_index_filename)) {
    Archive::Binary::In archive(contex_index_filename);
    try {
        RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

        archive.load(context_index, progress_bar, "context index");
    } catch (RACES::Archive::WrongFileFormatDescr& ex) {
        raise_error(ex, "context index");
    } catch (RACES::Archive::WrongFileFormatVersion& ex) {
        raise_error(ex, "context index");
    }

    return context_index;
  }

  using namespace Rcpp;

  Rcout << "Building context index..." << std::endl << std::flush;

  std::set<GenomicRegion> regions_to_avoid;

  auto drivers_path = storage.get_driver_mutations_path();
  if (std::filesystem::exists(drivers_path)) {
    auto driver_storage = DriverStorage::load(drivers_path);

    for (const auto& [name, mutation] : driver_storage.get_mutations()) {
        regions_to_avoid.emplace(mutation, std::max(static_cast<size_t>(1),
                                                    mutation.ref.size()));
    }
  }

  std::list<GenomicRegion> chr_regions;
  {
    UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);
    const auto reference_path = storage.get_reference_path();
    context_index = Index::build_index(reference_path, regions_to_avoid,
                                       context_sampling, &progress_bar);
    chr_regions = context_index.get_chromosome_regions();
  }

  if (chr_regions.size() > 0) {
    UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);
    Archive::Binary::Out archive(contex_index_filename);

    archive.save(context_index, progress_bar, "context index");
  }

  Rcout << "done" << std::endl;

  return context_index;
}

inline std::filesystem::path
get_rs_index_path(const GenomicDataStorage& storage,
                  const size_t max_motif_size,
                  const size_t max_repetition_storage)
{
    return storage.get_directory()/std::string("rs_index_"
                                                + std::to_string(max_motif_size)
                                                + "_"
                                                + std::to_string(max_repetition_storage)
                                                + ".rsif");
}

RACES::Mutations::RSIndex
build_rs_index(const GenomicDataStorage& storage,
               const size_t max_motif_size,
               const size_t max_repetition_storage,
               const bool& quiet)
{
  using namespace RACES;
  using namespace RACES::Mutations;

  using Index = RSIndex;

  Index rs_index;

  auto rs_index_filename = get_rs_index_path(storage, max_motif_size,
                                             max_repetition_storage);

  if (std::filesystem::exists(rs_index_filename)) {
    Archive::Binary::In archive(rs_index_filename);
    try {
        RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

        archive.load(rs_index, progress_bar, "RS index");
    } catch (RACES::Archive::WrongFileFormatDescr& ex) {
        raise_error(ex, "RS index");
    } catch (RACES::Archive::WrongFileFormatVersion& ex) {
        raise_error(ex, "RS index");
    }

    return rs_index;
  }

  using namespace Rcpp;

  if (!quiet) {
    Rcout << "Building repeated sequence index..."
          << std::endl << std::flush;
  }

/*
  std::set<GenomicRegion> regions_to_avoid;

  auto drivers_path = storage.get_driver_mutations_path();
  if (std::filesystem::exists(drivers_path)) {
    auto driver_storage = DriverStorage::load(drivers_path);

    for (const auto& [name, mutation] : driver_storage.get_mutations()) {
        regions_to_avoid.emplace(mutation, std::max(static_cast<size_t>(1),
                                                    mutation.ref.size()));
    }
  }
*/

  {
    UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);
    const auto reference_path = storage.get_reference_path();
    rs_index = Index::build_index(reference_path, max_motif_size,
                                  max_repetition_storage, &progress_bar);
  }

  Archive::Binary::Out archive(rs_index_filename);

  UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

  archive.save(rs_index, progress_bar, "RS index");

  if (!quiet) {
    Rcout << "done" << std::endl;
  }

  return rs_index;
}

template<typename ABSOLUTE_GENOTYPE_POSITION>
std::map<RACES::Mutations::ChromosomeId, size_t>
get_num_of_alleles(const RACES::Mutations::ContextIndex<ABSOLUTE_GENOTYPE_POSITION>& context_index,
                   const size_t& default_num_of_alleles,
                   const std::map<std::string, size_t>& alleles_num_exceptions)
{
  using namespace RACES::Mutations;

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


template<typename MUTATION_TYPE,
        std::enable_if_t<std::is_base_of_v<RACES::Mutations::MutationType, MUTATION_TYPE>, bool> = true>
std::map<std::string, RACES::Mutations::Signature<MUTATION_TYPE>>
load_signature(const GenomicDataStorage& storage)
{
  std::ifstream is(storage.get_signatures_path<MUTATION_TYPE>());

  return RACES::Mutations::Signature<MUTATION_TYPE>::read_from_stream(is);
}

RACES::Mutations::GenomicRegion get_CNA_region(const RACES::IO::CSVReader::CSVRow& row, const size_t& row_num)
{
  using namespace RACES::Mutations;

  ChromosomeId chr_id;
  try {
    chr_id = GenomicPosition::stochr(row.get_field(0));
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


std::vector<RACES::Mutations::CNA> load_passenger_CNAs(const std::filesystem::path& CNAs_csv,
                                                       const std::string& tumour_type,
                                                       const std::string& tumour_study)
{
  std::set<RACES::Mutations::CNA> CNAs;

  RACES::IO::CSVReader csv_reader(CNAs_csv, true, '\t');

  size_t row_num{2};
  for (const auto& row : csv_reader) {
    if (row.size()<7) {
      throw std::runtime_error("The CNA CSV must contains at least 7 columns");
    }
    if (((tumour_type=="") || (row.get_field(5) == tumour_type))
        && ((tumour_study=="") || (row.get_field(6) == tumour_study))) {
      const auto region = get_CNA_region(row, row_num);

      const auto major = row.get_field(3);
      try {
        if (major=="NA" || (stoi(major)>1)) {
          CNAs.emplace(region.get_begin(), region.size(),
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
          CNAs.emplace(region.get_begin(), region.size(),
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

  return std::vector<RACES::Mutations::CNA>(CNAs.begin(), CNAs.end());
}

void MutationEngine::init_mutation_engine(const bool& quiet)
{
  context_index = build_contex_index(storage, context_sampling, quiet);
  rs_index = build_rs_index(storage, max_motif_size,
                            max_repetition_storage, quiet);

  reset();
}

MutationEngine::MutationEngine(const std::string& setup_name,
                               const std::string& germline_subject,
                               const size_t& context_sampling,
                               const size_t& max_motif_size,
                               const size_t& max_repetition_storage,
                               const std::string& tumour_type,
                               const std::string& tumor_study,
                               const bool& quiet):
  storage(setup_storage(setup_name)), germline_subject(germline_subject),
  context_sampling(context_sampling), max_motif_size(max_motif_size),
  max_repetition_storage(max_repetition_storage), tumour_type(tumour_type),
  tumor_study(tumor_study)
{
  auto setup_cfg = supported_setups.at(setup_name);

  init_mutation_engine(quiet);
}

MutationEngine::MutationEngine(const std::string& directory,
                               const std::string& reference_source,
                               const std::string& SBS_signatures_source,
                               const std::string& indel_signatures_source,
                               const std::string& drivers_source,
                               const std::string& passenger_CNAs_source,
                               const std::string& germline_source,
                               const std::string& germline_subject,
                               const size_t& context_sampling,
                               const size_t& max_motif_size,
                               const size_t& max_repetition_storage,
                               const std::string& tumour_type,
                               const std::string& tumor_study,
                               const bool& quiet):
  storage(setup_storage(directory, reference_source, SBS_signatures_source,
                        indel_signatures_source, drivers_source,
                        passenger_CNAs_source, germline_source)),
  germline_subject(germline_subject), context_sampling(context_sampling),
  max_motif_size(max_motif_size),
  max_repetition_storage(max_repetition_storage),
  tumour_type(tumour_type), tumor_study(tumor_study)
{
  init_mutation_engine(quiet);
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
                                     const std::string& SBS_signatures_source,
                                     const std::string& indel_signatures_source,
                                     const std::string& drivers_source,
                                     const std::string& passenger_CNAs_source,
                                     const std::string& germline_source,
                                     const std::string& setup_code,
                                     const std::string& germline_subject,
                                     const size_t& context_sampling,
                                     const size_t& max_motif_size,
                                     const size_t& max_repetition_storage,
                                     const std::string& tumour_type,
                                     const std::string& tumor_study,
                                     const bool quiet)
{
  if (setup_code!="") {
    if (directory!="" || reference_source!="" || SBS_signatures_source!=""
         || drivers_source!="" || passenger_CNAs_source !=""
         || germline_source !="") {
      throw std::domain_error("when \"setup_code\" is NOT provided, the parameters "
                              "\"directory\", \"reference_src\", \"SBS_signatures_src\", "
                              "\"SBS_signatures_src\", \"passenger_CNAs_src\", and "
                              "\"germline_src\" must be avoided.");
    }

    return MutationEngine(setup_code, germline_subject, context_sampling,
                          max_motif_size, max_repetition_storage, tumour_type,
                          tumor_study, quiet);
  }

  if (directory=="" || reference_source=="" || SBS_signatures_source==""
      || indel_signatures_source== "" || passenger_CNAs_source== ""
      || germline_source== "") {
    throw std::domain_error("when \"setup_code\" is NOT provided, the parameters "
                            "\"directory\", \"reference_src\", \"SBS_signatures_src\", "
                            "\"SBS_signatures_src\", \"passenger_CNAs_src\", and "
                            "\"germline_src\" are mandatory.");
  }

  return MutationEngine(directory, reference_source, SBS_signatures_source,
                        indel_signatures_source, drivers_source,
                        passenger_CNAs_source, germline_source, germline_subject,
                        context_sampling, max_motif_size, max_repetition_storage,
                        tumour_type, tumor_study, quiet);

}

Rcpp::List MutationEngine::get_available_tumour_type(const std::string& setup_code)
{
    auto setup_cfg = supported_setups.at(setup_code);

    auto storage = setup_storage(setup_code);

    std::set<std::pair<std::string, std::string>> CNA_type_study;

    RACES::IO::CSVReader csv_reader(storage.get_passenger_CNAs_path(), true, '\t');
    for (const auto& row : csv_reader) {
        if (row.size()<7) {
            throw std::runtime_error("The CNA CSV must contains at least 7 columns");
        }

        CNA_type_study.insert({row.get_field(5), row.get_field(6)});
    }

    using namespace Rcpp;

    StringVector types(CNA_type_study.size()), studies(CNA_type_study.size());

    size_t i{0};
    for (const auto& [type, study] : CNA_type_study) {
        types[i] = type;
        studies[i] = study;

        ++i;
    }
    
    return DataFrame::create(_["type"]=types, _["study"]=studies);
}


void MutationEngine::add_exposure(const double& time, const Rcpp::List& exposure)
{
  auto c_exposure = get_map<double, TestNonNegative>(exposure);

  m_engine.add(time, c_exposure);
}

const RACES::Mutations::SID&
get_mutation_from_name(const std::map<std::string, RACES::Mutations::SID>& driver_code_map,
                       const Rcpp::CharacterVector& R_mutation_code)
{
    const auto mutation_code = Rcpp::as<std::string>(R_mutation_code);

    const auto found = driver_code_map.find(mutation_code);

    if (found == driver_code_map.end()) {
        throw std::domain_error("Unknown mutation code " + mutation_code + ".");
    }

    return found->second;
}

SIDSpec get_mutation_from_list(const std::map<std::string, RACES::Mutations::SID>& driver_code_map,
                               const Rcpp::List& SID_spec, const size_t index)
{
    const size_t spec_size = static_cast<size_t>(SID_spec.size());
    if (spec_size > 2 || spec_size == 0) {
        throw std::domain_error("The " + ordtostr(index)
                                + " element in the driver mutation list"
                                + " is not an mutation specification");
    }

    if (TYPEOF(SID_spec[0])!=STRSXP) {
        throw std::domain_error("The " + ordtostr(index)
                                + " element in the driver mutation list"
                                + " is not an mutation specification");
    }

    RACES::Mutations::AlleleId allele_id = RANDOM_ALLELE;

    if (spec_size>1) {
        if (TYPEOF(SID_spec[1])!=REALSXP) {
            throw std::domain_error("The " + ordtostr(index)
                                + " element in the driver mutation list"
                                + " is not an mutation specification");
        }

        allele_id = Rcpp::as<RACES::Mutations::AlleleId>(SID_spec[1]);
    }

    return SIDSpec(allele_id, get_mutation_from_name(driver_code_map, SID_spec[0]));
}

void
get_mutation_spec(std::list<SIDSpec>& c_sids,
                  std::list<RACES::Mutations::CNA>& c_cnas,
                  std::list<RACES::Mutations::DriverMutations::MutationType>& application_order,
                  const std::map<std::string, RACES::Mutations::SID>& driver_code_map,
                  const Rcpp::List& rcpp_list, const size_t& index)
{
    switch (TYPEOF(rcpp_list[index])) {
        case STRSXP:
            c_sids.emplace_back(RANDOM_ALLELE,
                                get_mutation_from_name(driver_code_map,
                                                       rcpp_list[index]));

            application_order.push_back(RACES::Mutations::DriverMutations::SID_TURN);

            return;
        case VECSXP:
            c_sids.push_back(get_mutation_from_list(driver_code_map,
                                                    rcpp_list[index],
                                                    index+1));

            application_order.push_back(RACES::Mutations::DriverMutations::SID_TURN);

            return;
        case S4SXP:
            {
                Rcpp::S4 s4obj( rcpp_list[index] );
                if ( s4obj.is("Rcpp_Mutation")) {
                    const auto sid = Rcpp::as<SID>(rcpp_list[index]);

                    c_sids.push_back(static_cast<SIDSpec>(sid));

                    application_order.push_back(RACES::Mutations::DriverMutations::SID_TURN);

                    return;
                }
                if ( s4obj.is("Rcpp_CNA")) {
                    const auto cna = Rcpp::as<CNA>(rcpp_list[index]);

                    c_cnas.push_back(static_cast<const RACES::Mutations::CNA&>(cna));

                    application_order.push_back(RACES::Mutations::DriverMutations::CNA_TURN);

                    return;
                }
                if ( s4obj.is("Rcpp_WholeGenomeDoubling")) {
                    application_order.push_back(RACES::Mutations::DriverMutations::WGD_TURN);

                    return;
                }
            }
            break;
        default:
            break;
    }
    std::cout << TYPEOF(rcpp_list[index]) << std::endl;
    throw std::domain_error("The " + ordtostr(index+1)
                            + " element in the driver mutation list"
                            + " is not an mutation specification");
}

void
get_mutation_lists(std::list<SIDSpec>& c_sids,
                   std::list<RACES::Mutations::CNA>& c_cnas,
                   std::list<RACES::Mutations::DriverMutations::MutationType>& application_order,
                   const std::map<std::string, RACES::Mutations::SID>& driver_code_map,
                   const Rcpp::List& rcpp_list)
{
    const size_t list_size = static_cast<size_t>(rcpp_list.size());
    for (size_t i=0; i<list_size; ++i) {
        get_mutation_spec(c_sids, c_cnas, application_order,
                          driver_code_map, rcpp_list, i);
    }
}

void MutationEngine::add_mutant(const std::string& mutant_name,
                                const Rcpp::List& epistate_passenger_rates)
{
  Rcpp::List empty_list;

  add_mutant(mutant_name, epistate_passenger_rates, empty_list);
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
    if (names[i]!="SNV" && names[i]!="CNA" && names[i]!="indel") {
      return false;
    }
  }

  return true;
}

RACES::Mutations::PassengerRates
get_passenger_rates(const Rcpp::List& list)
{
  RACES::Mutations::PassengerRates p_rates;

  if (!list.hasAttribute("names")) {
    throw std::runtime_error("Passenger rates list must be a named list whose names "
                             "are in the set {\"SNV\", \"indel\", \"CNA\"}.");
  }

  Rcpp::CharacterVector names = list.names();
  const size_t list_size = static_cast<size_t>(list.size());
  for (size_t i=0; i<list_size; ++i) {
    if (names[i]=="SNV") {
      p_rates.snv = get_non_negative(list[i], ": SNV rates must be non-negative");
    } else if (names[i]=="CNA") {
      p_rates.cna = get_non_negative(list[i], ": CNA rates must be non-negative");
    } else if (names[i]=="indel") {
      p_rates.indel = get_non_negative(list[i], ": indel rates must be non-negative");
    } else {
        throw std::runtime_error("\"" + names[i] + "\" is an unsupported name in "
                                 + "passenger rates list.");
    }
  }

  return p_rates;
}

std::map<std::string, RACES::Mutations::PassengerRates>
get_epistate_passenger_rates(const Rcpp::List& list)
{
  std::map<std::string, RACES::Mutations::PassengerRates> ep_rates;

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

struct FilterNonChromosomeSequence : public RACES::IO::FASTA::SequenceFilter
{
    RACES::Mutations::ChromosomeId last_chr_id;

    inline bool operator()(const std::string& header)
    {
        return !RACES::IO::FASTA::is_chromosome_header(header, last_chr_id);
    }
};

void check_wrong_chromosome_SNV(const std::map<RACES::Mutations::ChromosomeId, SID_iterator>& SNV_partition)
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
                                 const std::filesystem::path fasta_filename,
                                 std::list<SIDSpec>& SNVs)
{
  RACES::UI::ProgressBar progress_bar(Rcpp::Rcout);

  std::map<RACES::Mutations::ChromosomeId, SID_iterator> SNV_partition;

  size_t SNV_to_check{0};
  for (auto it=SNVs.begin(); it != SNVs.end(); ++it) {
    if (it->ref == "?") {
      SNV_partition[it->chr_id].push_back(it);
      ++SNV_to_check;
    }
  }

  RACES::IO::FASTA::Sequence chr_seq;

  FilterNonChromosomeSequence filter;

  progress_bar.set_message("Retrieving \""+ mutant_name + "\" SNVs");

  const auto fasta_size = filesize(fasta_filename);
  std::ifstream fasta_stream(fasta_filename);
  while (SNV_to_check>0 && RACES::IO::FASTA::Sequence::read(fasta_stream, chr_seq, filter, progress_bar)) {
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

        if (SNV_it->ref == "?") {
          const auto& candidate_ref = chr_seq.nucleotides[SNV_it->position-1];

          SNV_it->ref = std::string(1, candidate_ref);
        }
        --SNV_to_check;
      }

      SNV_partition.erase(found);
    }
  }

  check_wrong_chromosome_SNV(SNV_partition);

  progress_bar.set_progress(100, "\"" + mutant_name + "\" SNVs retrieved");
}

#include <iostream>

void MutationEngine::add_mutant(const std::string& mutant_name,
                                const Rcpp::List& epistate_passenger_rates,
                                const Rcpp::List& drivers)
{
  std::list<SIDSpec> c_sids;
  std::list<RACES::Mutations::CNA> c_cnas;

  const auto& driver_storage = m_engine.get_driver_storage();

  std::list<RACES::Mutations::DriverMutations::MutationType> application_order;

  get_mutation_lists(c_sids, c_cnas, application_order,
                     driver_storage.get_mutations(), drivers);

  retrieve_missing_references(mutant_name, storage.get_reference_path(), c_sids);

  if (contains_passenger_rates(epistate_passenger_rates)) {
    auto p_rates = get_passenger_rates(epistate_passenger_rates);

    m_engine.add_mutant(mutant_name, {{"", p_rates}}, c_sids, c_cnas,
                        application_order);

    return;
  }

  auto epi_rates = get_epistate_passenger_rates(epistate_passenger_rates);
  m_engine.add_mutant(mutant_name, epi_rates, c_sids, c_cnas,
                      application_order);
}

PhylogeneticForest
MutationEngine::place_mutations(const SamplesForest& forest,
                                const size_t& num_of_preneoplatic_SNVs,
                                const std::string& preneoplatic_SNV_signature_name,
                                const size_t& num_of_preneoplatic_indels,
                                const std::string& preneoplatic_indel_signature_name,
                                const int seed)
{
  RACES::UI::ProgressBar progress_bar(Rcpp::Rcout);

  progress_bar.set_message("Placing mutations");

  auto phylo_forest = m_engine.place_mutations(forest, num_of_preneoplatic_SNVs,
                                               num_of_preneoplatic_indels,
                                               progress_bar, seed,
                                               preneoplatic_SNV_signature_name,
                                               preneoplatic_indel_signature_name);
  progress_bar.set_message("Mutations placed");

  using MutationType = RACES::Mutations::MutationType;

  const auto& const_m_engine = m_engine;

  return {std::move(phylo_forest), germline_subject, storage.get_reference_path(),
          const_m_engine.get_timed_exposures(MutationType::Type::SBS),
          const_m_engine.get_timed_exposures(MutationType::Type::INDEL)};
}

Rcpp::List MutationEngine::get_SNV_signatures_dataframe() const
{
  Rcpp::Function read_delim("read.delim");

  return read_delim(to_string(storage.get_signatures_path<RACES::Mutations::SBSType>()),
                    Rcpp::_["quote"]="");
}

Rcpp::List MutationEngine::get_indel_signatures_dataframe() const
{
  Rcpp::Function read_delim("read.delim");

  return read_delim(to_string(storage.get_signatures_path<RACES::Mutations::IDType>()),
                    Rcpp::_["quote"]="");
}

Rcpp::List MutationEngine::get_known_driver_mutations() const
{
  auto driver_filename = m_engine.get_driver_storage().get_source_path().string();

  Rcpp::Function read_delim("read.delim");

  return read_delim(driver_filename, Rcpp::_["quote"]="");
}

template<typename ITERATOR>
std::ostream& show_list(std::ostream& os, ITERATOR it, ITERATOR last, const std::string& front="")
{
  for (; it != last; ++it) {
    os << front << *it << std::endl;
  }

  return os;
}

std::ostream&
show_driver_mutations(std::ostream& os,
                      const RACES::Mutations::DriverMutations& driver_mutations, 
                      const std::string& indent="")
{
    using namespace RACES::Mutations;

    auto SID_it = driver_mutations.SIDs.begin();
    auto CNA_it = driver_mutations.CNAs.begin();

    for (auto order_it = driver_mutations.application_order.begin();
        order_it != driver_mutations.application_order.end(); ++order_it) {

        switch(*order_it) {
            case DriverMutations::SID_TURN:
                os << indent << *SID_it << std::endl;
                ++SID_it;

                break;
            case DriverMutations::CNA_TURN:
                os << indent << *CNA_it << std::endl;
                ++CNA_it;

                break;
            case DriverMutations::WGD_TURN:
                os << indent << "Whole genome duplication" << std::endl;
                break;
            default:
                throw std::runtime_error("Unsupported driver mutation type.");
        }
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
    if (p_rates.indel>0) {
      Rcout << sep << "indel: " << p_rates.indel;
      sep = ", ";
    }
    if (p_rates.cna>0) {
      Rcout << sep << "CNA: " << p_rates.cna;
    }
    Rcout << "}";
  }

  Rcout << std::endl << std::endl << " Driver mutations" << std::endl;
  for (const auto&[mutant_name, driver_mutations]: m_properties.get_driver_mutations()) {
    Rcout << "   \"" << mutant_name << "\":" << std::endl;
    show_driver_mutations(Rcout, driver_mutations, "       ");
    if (driver_mutations.application_order.size()==0) {
      Rcout << "   No driver mutations for \"" << mutant_name << "\"" << std::endl;
    }
  }

  Rcout << std::endl << " Timed Exposure" << std::endl;

  show_timed_exposures<RACES::Mutations::SBSType>();
  Rcout << std::endl;
  show_timed_exposures<RACES::Mutations::IDType>();
  Rcout << std::endl;
}

void MutationEngine::rebuild_indices(const bool quiet)
{
  auto index_path = get_context_index_path(storage, context_sampling);

  std::filesystem::remove(index_path);

  context_index = build_contex_index(storage, quiet, context_sampling);

  index_path = get_rs_index_path(storage, max_motif_size,
                                 max_repetition_storage);

  std::filesystem::remove(index_path);

  rs_index = build_rs_index(storage, max_motif_size,
                            max_repetition_storage, quiet);
}

void MutationEngine::set_context_sampling(const size_t& context_sampling,
                                          const bool quiet)
{
  this->context_sampling = context_sampling;

  context_index = build_contex_index(storage, context_sampling, quiet);

  reset();
}

void MutationEngine::reset(const bool full, const bool quiet)
{
  using namespace RACES::Mutations;

  MutationalProperties mutational_properties;
  std::map<MutationType::Type, std::map<RACES::Time, MutationalExposure>> timed_exposures;

  if (!full) {
    mutational_properties = m_engine.get_mutational_properties();

    const auto& const_m_engine = m_engine;

    timed_exposures[MutationType::Type::SBS] = const_m_engine.get_timed_exposures(MutationType::Type::SBS);
    timed_exposures[MutationType::Type::INDEL] = const_m_engine.get_timed_exposures(MutationType::Type::INDEL);
  }

  auto SBS_signatures = load_signature<SBSType>(storage);
  auto indel_signatures = load_signature<IDType>(storage);

  auto passenger_CNAs = load_passenger_CNAs(storage.get_passenger_CNAs_path(),
                                            tumour_type, tumor_study);

  auto driver_storage = DriverStorage::load(storage.get_driver_mutations_path());

  auto& germline_storage = storage.get_germline_storage();

  if (germline_subject == "") {
    auto germline_subjects = germline_storage.get_population();

    if (germline_subjects.size() == 0) {
      throw std::runtime_error("No germline subject available.");
    }

    germline_subject = germline_subjects[0].name;
  }

  auto germline = germline_storage.get_germline(germline_subject, quiet);

  m_engine = RACES::Mutations::MutationEngine(context_index, rs_index,
                                              SBS_signatures,
                                              indel_signatures,
                                              mutational_properties, germline,
                                              driver_storage, passenger_CNAs);

  for (const auto& [type, mutation_timed_exposures] : timed_exposures) {
    for (const auto& [time, exposure] : mutation_timed_exposures) {
      m_engine.add(time, exposure);
    }
  }
}

void MutationEngine::set_germline_subject(const std::string& germline_subject)
{
  storage.get_germline_storage().get_subject(germline_subject);

  this->germline_subject = germline_subject;

  reset(false);
}

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

#include <Rcpp.h>


#include <germline.hpp>

#include <csv_reader.hpp>
#include <context_index.hpp>
#include <progress_bar.hpp>
#include <utils.hpp>


#include "genomic_data_storage.hpp"


GermlineSubject::GermlineSubject(const std::string& name, const std::string& population,
                                 const std::string& super_population, const std::string& gender):
  name(name), population(population), super_population(super_population), gender(gender)
{}

GermlineStorage::GermlineStorage()
{}

GermlineStorage::GermlineStorage(const std::filesystem::path& directory):
  directory(directory)
{
  if (!std::filesystem::exists(directory)) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not exist.");
  }
  if (!std::filesystem::is_directory(directory)) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" is not a "
                             + "directory.");
  }
  if (!std::filesystem::exists(get_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"germline.csv\".");
  }
  if (!std::filesystem::exists(get_population_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"population.csv\".");
  }
  if (!std::filesystem::exists(get_population_descriptions_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"population_descr"
                             + "iptions.csv\".");
  }
  if (!std::filesystem::exists(get_alleles_file())) {
    throw std::runtime_error("Designed germline mutations directory \""
                             + to_string(directory) + "\" does not "
                             + "contains the file \"alleles_per_chr."
                             + "csv\".");
  }
}

std::vector<GermlineSubject> GermlineStorage::get_population() const
{
  std::list<GermlineSubject> s_list;

  Races::IO::CSVReader csv_reader(get_population_file(), true, '\t');

  for (const auto& row : csv_reader) {
    s_list.emplace_back(row.get_field(0), row.get_field(1),
                        row.get_field(2), row.get_field(3));
  }

  return {s_list.begin(), s_list.end()};
}

Rcpp::List GermlineStorage::get_population_df() const
{
  using namespace Rcpp;

  Function read_csv("read.csv");

  return read_csv(_["file"]=to_string(get_population_file()),
                  _["quote"]="", _["header"]=true, _["sep"] = "\t");
}

Rcpp::List GermlineStorage::get_population_descritions_df() const
{
  using namespace Rcpp;

  Function read_csv("read.csv");

  return read_csv(_["file"]=get_population_descriptions_file(),
                  _["quote"]="", _["header"]=true, _["sep"] = "\t");
}

std::map<Races::Mutations::ChromosomeId, size_t>
GermlineStorage::get_alleles_per_chromosome(const std::string& gender) const
{
  using namespace Races::Mutations;

  std::map<ChromosomeId, size_t> alleles_per_chromosome;

  Races::IO::CSVReader csv_reader(get_alleles_file(), true, '\t');

  const auto& header = csv_reader.get_header();
  auto found = find(header.begin(), header.end(), gender);

  if (found == header.end()) {
    throw std::runtime_error("Unknown gender " + gender + ".");
  }

  size_t index = found-header.begin();
  for (const auto& row : csv_reader) {
    auto chr_id = GenomicPosition::stochr(row.get_field(0));
    alleles_per_chromosome.insert({chr_id,
                                   std::stoul(row.get_field(index))});
  }

  return alleles_per_chromosome;
}

GermlineSubject GermlineStorage::get_subject(const std::string& subject_name) const
{
  auto gemline_subjects = get_population();

  if (gemline_subjects.size()==0) {
    throw std::runtime_error("No germline subject available.");
  }

  for (const auto& subject : gemline_subjects) {
    if (subject_name == subject.name) {
      return subject;
    }
  }

  throw std::runtime_error("Germline subject \"" + subject_name
                           + "\" not available.");
}


Rcpp::List GermlineStorage::get_subject_df(const std::string& subject_name) const
{
  using namespace Rcpp;

  Races::IO::CSVReader csv_reader(get_population_file(), true, '\t');

  for (const auto& row : csv_reader) {
    if (row.get_field(0) == subject_name) {
      CharacterVector sample(1), pop(1), super_pop(1), gender(1);

      sample[0] = row.get_field(0);
      pop[0] = row.get_field(1);
      super_pop[0] = row.get_field(2);
      gender[0] = row.get_field(3);

      return DataFrame::create(_["sample"] = sample, _["pop"] = pop,
                               _["super_pop"] = super_pop,
                               _["gender"] = gender);
    }
  }

  throw std::runtime_error("Germline subject \"" + subject_name
                           + "\" unknown.");
}

Races::Mutations::GenomeMutations
GermlineStorage::build_germline(const std::string& subject_name) const
{
  using namespace Races::Mutations;

  auto bin_path = get_binary_file(subject_name);

  auto subject = get_subject(subject_name);

  auto num_of_alleles = get_alleles_per_chromosome(subject.gender);

  auto germline = GermlineMutations::load(get_file(), num_of_alleles,
                                          subject.name);

  Races::Archive::Binary::Out oarchive(bin_path);

  oarchive.save(germline, "germline");

  return germline;
}

Races::Mutations::GenomeMutations
GermlineStorage::get_germline(const std::string& subject_name) const
{
  using namespace Races::Mutations;

  auto bin_path = get_binary_file(subject_name);

  if (!std::filesystem::exists(bin_path)) {
    return build_germline(subject_name);
  }

  Races::Archive::Binary::In iarchive(bin_path);

  GenomeMutations germline;

  iarchive.load(germline, "germline");

  return germline;
}

bool is_an_URL(const std::string& reference)
{
  std::set<std::string> protocols{"ftp", "http"};

  for (const auto& protocol : protocols) {
    if (reference.find(protocol)==0) {
      return true;
    }
  }

  return false;
}

GenomicDataStorage::GenomicDataStorage(const std::string& directory,
                                       const std::string& reference_source,
                                       const std::string& SBS_source,
                                       const std::string& driver_mutations_source,
                                       const std::string& passenger_CNAs_source,
                                       const std::string& germline_source):
  directory(directory), reference_src(reference_source),
  SBS_src(SBS_source), drivers_src(driver_mutations_source),
  passenger_CNAs_src(passenger_CNAs_source), germline_src(germline_source)
{
  std::filesystem::create_directory(directory);

  retrieve_reference();
  retrieve_SBS();
  retrieve_drivers();
  retrieve_passenger_CNAs();
  const auto germline_path = retrieve_germline();

  germline_storage = GermlineStorage(germline_path);
}

std::filesystem::path GenomicDataStorage::get_destination_path(const std::string& url) const
{
  try {
    auto name_occurrence = url.find_last_of('/')+1;

    auto filename = url.substr(name_occurrence);

    auto question_occurrence = filename.find('?');

    if (question_occurrence != std::string::npos) {
      filename = filename.substr(0, question_occurrence);
    }

    return directory/filename;

  } catch (std::exception& e) {
    throw std::domain_error("\""+url+"\" is not a valid URL.");
  }
}

std::filesystem::path GenomicDataStorage::download_file(const std::string& url) const
{
  if (!std::filesystem::exists(directory)) {
    std::filesystem::create_directory(directory);
  }

  using namespace Rcpp;

  auto dest_filename = to_string(get_destination_path(url));

  // get default timeout option
  Function getOption_f("getOption");
  auto timeout = as<int>(getOption_f("timeout"));

  // raise the timeout to 1000 at least
  Function options_f("options");
  options_f(_["timeout"] = std::max(1000, timeout));

  // download the file
  Function download_f("download.file");
  download_f(_["url"] = url, _["destfile"] = dest_filename,
             _["mode"]="wb");

  // revert to the default timeout
  options_f(_["timeout"] = timeout);

  return dest_filename;
}

std::map<std::string, std::string> decompressors{
  {"gz", "gunzip"},
  {"bz2", "bunzip2"}
};

std::filesystem::path GenomicDataStorage::retrieve_reference()
{
  reference_downloaded = is_an_URL(reference_src);
  if (!reference_downloaded && !std::filesystem::exists(reference_src)) {
    throw std::runtime_error("Designed reference genome file \""
                        + reference_src + "\" does not exists.");
  }

  const auto reference_filename = get_reference_path();

  if (std::filesystem::exists(reference_filename)) {
    return reference_filename;
  }

  using namespace Rcpp;

  Rcout << "Downloading reference genome..." << std::endl << std::flush;

  auto downloaded_file = to_string(download_file(reference_src));

  Rcout << "Reference genome downloaded" << std::endl;

  auto suffix = downloaded_file.substr(downloaded_file.find_last_of('.')+1);

  if (suffix == "fa" && suffix == "fasta") {
    std::filesystem::rename(downloaded_file, reference_filename);
  } else {
    Rcout << "Decompressing reference file...";
    auto decomp_found = decompressors.find(suffix);

    if (decomp_found == decompressors.end()) {
      throw std::runtime_error("Unknown suffix \""+suffix+"\"");
    }

    Environment pkg = Environment::namespace_env("R.utils");
    Function decompress_f = pkg[decomp_found->second];

    decompress_f(_["filename"] = downloaded_file, 
                 _["destname"] = to_string(reference_filename));

    Rcout << "done" << std::endl;
  }

  return reference_filename;
}

std::filesystem::path GenomicDataStorage::retrieve_SBS()
{
  SBS_downloaded = is_an_URL(SBS_src);
  if (!SBS_downloaded && !std::filesystem::exists(SBS_src)) {
    throw std::runtime_error("Designed SBS file \"" + SBS_src
                             + "\" does not exists.");
  }

  const auto SBS_filename = get_SBS_path();

  if (std::filesystem::exists(SBS_filename)) {
    return SBS_filename;
  }

  using namespace Rcpp;

  Rcout << "Downloading SBS file..." << std::endl << std::flush;

  auto downloaded_file = download_file(SBS_src);

  Rcout << "SBS file downloaded" << std::endl;

  std::filesystem::rename(downloaded_file, SBS_filename);

  return SBS_filename;
}

std::filesystem::path GenomicDataStorage::retrieve_drivers()
{
  if (drivers_src == "") {
    return drivers_src;
  }

  drivers_downloaded = is_an_URL(drivers_src);
  if (!drivers_downloaded && !std::filesystem::exists(drivers_src)) {
    throw std::runtime_error("Designed driver mutations file \"" 
                             + drivers_src + "\" does not exists.");
  }

  const auto mutations_filename = get_driver_mutations_path();

  if (std::filesystem::exists(mutations_filename)) {
    return mutations_filename;
  }

  using namespace Rcpp;

  Rcout << "Downloading driver mutation file..." << std::endl << std::flush;

  auto downloaded_file = download_file(drivers_src);

  Rcout << "Driver mutation file downloaded" << std::endl;

  std::filesystem::rename(downloaded_file, mutations_filename);

  return mutations_filename;
}

std::filesystem::path GenomicDataStorage::retrieve_passenger_CNAs()
{
  passenger_CNAs_downloaded = is_an_URL(passenger_CNAs_src);
  if (!passenger_CNAs_downloaded 
       && !std::filesystem::exists(passenger_CNAs_src)) {
    throw std::runtime_error("Designed passenger CNAs file \"" + SBS_src
                             + "\" does not exists.");
  }

  const auto passenger_CNAs_filename = get_passenger_CNAs_path();

  if (std::filesystem::exists(passenger_CNAs_filename)) {
    return passenger_CNAs_filename;
  }

  using namespace Rcpp;

  Rcout << "Downloading passenger CNAs file..." << std::endl << std::flush;

  auto downloaded_file = download_file(passenger_CNAs_src);

  Rcout << "Passenger CNAs file downloaded" << std::endl;

  std::filesystem::rename(downloaded_file, passenger_CNAs_filename);

  return passenger_CNAs_filename;
}

std::filesystem::path GenomicDataStorage::retrieve_germline()
{
  germline_downloaded = is_an_URL(germline_src);
  if (!germline_downloaded) {
    if (!std::filesystem::exists(germline_src)) {
      throw std::runtime_error("Designed germline mutations directory \""
                               + germline_src + "\" does not exists.");
    }
    return germline_src;
  }

  const auto germline_path = directory/"germline_data";

  if (std::filesystem::exists(germline_path/"germlines.csv")) {
    return germline_path;
  }

  using namespace Rcpp;

  Rcout << "Downloading germline mutations..." << std::endl << std::flush;

  auto downloaded_file = download_file(germline_src);

  Rcout << "Germline mutations downloaded" << std::endl;

  Function untar("untar");
  untar(_["tarfile"] = to_string(downloaded_file),
        _["exdir"] = to_string(directory));

  return germline_path;
}

std::filesystem::path GenomicDataStorage::get_reference_path() const
{
  if (reference_downloaded) {
    return directory/std::string("reference.fasta");
  } else {
    return reference_src;
  }
}

std::filesystem::path GenomicDataStorage::get_SBS_path() const
{
  if (SBS_downloaded) {
    return directory/std::string("SBS.txt");
  } else {
    return SBS_src;
  }
}

std::filesystem::path GenomicDataStorage::get_passenger_CNAs_path() const
{
  if (passenger_CNAs_downloaded) {
    return directory/std::string("passenger_CNAs.txt");
  } else {
    return passenger_CNAs_src;
  }
}

std::filesystem::path GenomicDataStorage::get_driver_mutations_path() const
{
  if (drivers_downloaded) {
    return directory/std::string("drivers.txt");
  } else {
    return drivers_src;
  }
}

void GenomicDataStorage::save_sources() const
{
  std::filesystem::path param_file(directory/"sources.csv");
 
  std::ofstream of(directory/"sources.csv");

  of << "reference\t" << reference_src << std::endl
     << "SBS\t" << SBS_src << std::endl
     << "drivers\t" << drivers_src << std::endl
     << "passenger_CNAs\t" << passenger_CNAs_src << std::endl
     << "germline\t" << germline_src << std::endl;
}
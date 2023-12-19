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

#include "genomic_data_storage.hpp"


GenomicDataStorage::GenomicDataStorage(const std::string& directory):
  directory(directory)
{}

std::string GenomicDataStorage::get_destination_path(const std::string& url) const
{
  try {
    auto name_occurrence = url.find_last_of('/')+1;

    const auto filename = url.substr(name_occurrence);

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

  std::string dest_filename = get_destination_path(url);

  // get default timeout option
  Function getOption_f("getOption");
  auto timeout = as<int>(getOption_f("timeout"));

  // raise the timeout to 1000 at least 
  Function options_f("options");
  options_f(_["timeout"] = std::max(1000, timeout));

  // download the file
  Function download_f("download.file");
  download_f(_["url"] = url, _["destfile"] = dest_filename);

  // revert to the default timeout
  options_f(_["timeout"] = timeout);

  return dest_filename;
}

std::map<std::string, std::string> decompressors{
  {"gz", "gunzip"},
  {"bz2", "bunzip2"}
};

void GenomicDataStorage::download_reference(const std::string& url)
{
  std::string reference_filename = get_reference_path();

  if (std::filesystem::exists(reference_filename)) {
    return;
  }

  using namespace Rcpp;

  Rcout << "Downloading reference genome..." << std::endl << std::flush;

  std::string downloaded_file = download_file(url);

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

    decompress_f(_["filename"] = downloaded_file, _["destname"] = reference_filename);

    Rcout << "done" << std::endl;
  }

  reference_url = url;
}

void GenomicDataStorage::download_SBS(const std::string& url) 
{
  std::string SBS_filename = get_SBS_path();

  if (std::filesystem::exists(SBS_filename)) {
    return;
  }

  using namespace Rcpp;

  Rcout << "Downloading SBS file..." << std::endl << std::flush;

  std::string downloaded_file = download_file(url);

  Rcout << "SBS file downloaded" << std::endl;

  std::filesystem::rename(downloaded_file, SBS_filename);

  SBS_url = url;
}

void GenomicDataStorage::save_parameters() const
{
  std::filesystem::path param_file(directory/"parameters.txt");

  if (std::filesystem::exists(param_file)) {
    return;
  }
    
  std::ofstream of(directory/"parameters.txt");

  of << "reference_url: \"" << reference_url << "\"" << std::endl
     << "SBS_url: \"" << SBS_url << "\"" << std::endl;
}
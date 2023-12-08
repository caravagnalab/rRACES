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

#include "setup_mutation_engine.hpp"


class MutationGeneratorArranger
{
  std::filesystem::path directory;
  std::string reference_url;
  std::string SBS_url;
  size_t context_sampling;

  std::string get_destination_path(const std::string& url) const;

  std::filesystem::path download_file(const std::string& url) const;
public:
  MutationGeneratorArranger(const std::string& directory);

  void build_contex_index(const size_t context_sampling=100);

  void download_reference(const std::string& url);

  void download_SBS(const std::string& url);

  inline std::filesystem::path get_directory() const
  {
    return directory;
  }

  inline std::filesystem::path get_reference_path() const
  {
    return directory/std::string("reference.fasta");
  }

  inline std::filesystem::path get_SBS_path() const
  {
    return directory/std::string("SBS.txt");
  }

  inline std::filesystem::path get_context_index_path() const
  {
    return directory/std::string("context_index.cif");
  }

  void save_parameters() const;
};

MutationGeneratorArranger::MutationGeneratorArranger(const std::string& directory):
  directory(directory), context_sampling(0)
{}

std::string MutationGeneratorArranger::get_destination_path(const std::string& url) const
{
  try {
    auto name_occurrence = url.find_last_of('/')+1;

    const auto filename = url.substr(name_occurrence);

    return directory/filename;

  } catch (std::exception& e) {
    throw std::domain_error("\""+url+"\" is not a valid URL.");
  }  
}

std::filesystem::path MutationGeneratorArranger::download_file(const std::string& url) const
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

void MutationGeneratorArranger::download_reference(const std::string& url)
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


void MutationGeneratorArranger::download_SBS(const std::string& url) 
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

void MutationGeneratorArranger::build_contex_index(const size_t context_sampling)
{
  using namespace Races;
  using namespace Races::Passengers;

  std::string contex_index_filename = get_context_index_path();

  if (std::filesystem::exists(contex_index_filename)) {
    
    return;
  }

  using namespace Rcpp;

  Rcout << "Building context index..." << std::endl << std::flush;

  std::vector<GenomicRegion> chr_regions;
  {
    using Index = ContextIndex<uint32_t>;

    Index context_index;

    {
      UI::ProgressBar progress_bar;

      std::string reference_filename = get_reference_path();

      context_index = Index::build_index(reference_filename, context_sampling, &progress_bar);
      
      chr_regions = context_index.get_chromosome_regions();
    }

    if (chr_regions.size() > 0) {
      Archive::Binary::Out archive(contex_index_filename);

      archive.save(context_index, "context index");
    }

    Rcout << " Cleaning memory..." << std::flush;
  }

  Rcout << "done" << std::endl
        << "Context index built" << std::endl;
  this->context_sampling = context_sampling;
}

void MutationGeneratorArranger::save_parameters() const
{
  std::ofstream of(directory/"parameters.txt");

  of << "reference_url: \"" << reference_url << "\"" << std::endl
     << "SBS_url: \"" << SBS_url << "\"" << std::endl
     << "context_sampling: " << context_sampling << std::endl;
}

struct ArrangerSetting
{
  std::string description;
  std::filesystem::path directory;
  std::string reference_url;
  std::string SBS_url;

  ArrangerSetting(const std::string& description,
                  const std::filesystem::path& directory,
                  const std::string& reference_url,
                  const std::string& SBS_url):
    description(description), directory(directory),
    reference_url(reference_url), SBS_url(SBS_url)
  {}
};

std::map<std::string, ArrangerSetting> supported_codes{
  {
    "demo",
    {
      "A demostative set-up", "demo",
      "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2124/COSMIC_v3.4_SBS_GRCh38.txt"
    }
  },
  {
    "GRCh38",
    {
      "Homo sapiens (GRCh38)", "GRCh38",
      "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2124/COSMIC_v3.4_SBS_GRCh38.txt"
    }
  },
  {
    "GRCh37",
    {
      "Homo sapiens (GRCh37)", "GRCh37",
      "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
      "https://cancer.sanger.ac.uk/signatures/documents/2123/COSMIC_v3.4_SBS_GRCh37.txt"
    }
  }
};

void setup_mutation_engine_by_code(const std::string& setup_code,
                                   const size_t& context_sampling)
{
  using namespace Rcpp;

  auto code_it = supported_codes.find(setup_code);
  if (code_it == supported_codes.end()) {

    Rcerr << "\""+setup_code+"\" is an unknown code. "
          << "Supported codes are:" << std::endl;

    for (const auto& [s_code, setting]: supported_codes) {
      Rcerr << " - \""+s_code+"\": " << setting.description << std::endl;
    }

    return;
  }

  setup_mutation_engine(code_it->second.directory, 
                        code_it->second.reference_url,
                        code_it->second.SBS_url, context_sampling);
}

void setup_mutation_engine_by_code(const std::string& setup_code)
{
  setup_mutation_engine_by_code(setup_code, 100);
}

void setup_mutation_engine(const std::string& directory,
                           const std::string& reference_url,
                           const std::string& SBS_url,
                           const size_t& context_sampling)
{
  MutationGeneratorArranger arranger(directory);

  arranger.download_reference(reference_url);
  arranger.download_SBS(SBS_url);
  arranger.build_contex_index(context_sampling);

  arranger.save_parameters();
}

Rcpp::List get_mutation_engine_supported_codes()
{
  using namespace Rcpp;

  CharacterVector code_names(supported_codes.size()), code_descrs(supported_codes.size());

  size_t i=0;
  for (const auto& [s_code, setting]: supported_codes) {
    code_names[i] = s_code;
    code_descrs[i] = setting.description;
    ++i;
  }

  return DataFrame::create(_["code"]=code_names, _["description"]=code_descrs);
}

void setup_mutation_engine(const std::string& directory,
                           const std::string& reference_url,
                           const std::string& SBS_url)
{
  setup_mutation_engine(directory, reference_url, SBS_url, 100);
}

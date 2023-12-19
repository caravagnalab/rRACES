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

#ifndef __RRACES_SETUP_MUTATION_ENGINE__
#define __RRACES_SETUP_MUTATION_ENGINE__

#include <string>

#include <Rcpp.h>

#include <mutation_engine.hpp>


class GenomicDataStorage
{
  std::filesystem::path directory;
  std::string reference_url;
  std::string SBS_url;

  std::string get_destination_path(const std::string& url) const;

  std::filesystem::path download_file(const std::string& url) const;
public:
  GenomicDataStorage(const std::string& directory);

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

  void save_parameters() const;
};

#endif // __RRACES_SETUP_MUTATION_ENGINE__

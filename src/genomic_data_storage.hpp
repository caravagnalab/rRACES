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

#ifndef __RRACES_SETUP_MUTATION_ENGINE__
#define __RRACES_SETUP_MUTATION_ENGINE__

#include <string>

#include <Rcpp.h>

#include <mutation_engine.hpp>

struct GermlineSubject
{
  std::string name;
  std::string population;
  std::string super_population;
  std::string gender;

  GermlineSubject(const std::string& name, const std::string& population,
                  const std::string& super_population, const std::string& gender);
};

class GermlineStorage
{
  std::filesystem::path directory;

  inline std::filesystem::path get_alleles_file() const
  {
    return get_path()/std::string("alleles_per_chr.csv");
  }

  inline std::filesystem::path get_population_file() const
  {
    return get_path()/std::string("population.csv");
  }

  inline std::filesystem::path get_population_descriptions_file() const
  {
    return get_path()/std::string("population_descriptions.csv");
  }

  inline std::filesystem::path get_file() const
  {
    return get_path()/std::string("germlines.csv");
  }

  inline std::filesystem::path get_binary_file(const std::string& subject_name) const
  {
    return get_path()/("germline_" + subject_name + ".dat");
  }

  Races::Mutations::GenomeMutations build_germline(const std::string& subject_name) const;

public:

  GermlineStorage();

  GermlineStorage(const std::filesystem::path& directory);

  inline std::filesystem::path get_path() const
  {
    return directory;
  }

  std::vector<GermlineSubject> get_population() const;

  std::map<Races::Mutations::ChromosomeId, size_t>
  get_alleles_per_chromosome(const std::string& gender) const;

  GermlineSubject get_subject(const std::string& subject_name) const;

  Races::Mutations::GenomeMutations get_germline(const std::string& subject_name) const;

  Rcpp::List get_subject_df(const std::string& subject_name) const;
  
  Rcpp::List get_population_df() const;
  
  Rcpp::List get_population_descritions_df() const;
};

class GenomicDataStorage
{
  std::filesystem::path directory;
  GermlineStorage germline_storage;

  bool reference_downloaded;
  std::string reference_src;

  bool SBS_downloaded;
  std::string SBS_src;

  bool drivers_downloaded;
  std::string drivers_src;

  bool passenger_CNAs_downloaded;
  std::string passenger_CNAs_src;

  bool germline_downloaded;
  std::string germline_src;

  std::string get_destination_path(const std::string& url) const;

  std::filesystem::path download_file(const std::string& url) const;

  std::filesystem::path retrieve_reference();

  std::filesystem::path retrieve_SBS();

  std::filesystem::path retrieve_drivers();

  std::filesystem::path retrieve_passenger_CNAs();

  std::filesystem::path retrieve_germline();
public:
  GenomicDataStorage(const std::string& directory,
                     const std::string& reference_source,
                     const std::string& SBS_source,
                     const std::string& driver_mutations_source,
                     const std::string& passenger_CNAs_source,
                     const std::string& germline_source);

  inline std::filesystem::path get_directory() const
  {
    return directory;
  }

  std::filesystem::path get_reference_path() const;

  std::filesystem::path get_SBS_path() const;

  std::filesystem::path get_driver_mutations_path() const;

  std::filesystem::path get_passenger_CNAs_path() const;

  inline const GermlineStorage& get_germline_storage() const
  {
    return germline_storage;
  }

  void save_sources() const;
};

#endif // __RRACES_SETUP_MUTATION_ENGINE__

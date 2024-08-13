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

#ifndef __RRACES_PHYLOGENETIC_FOREST__
#define __RRACES_PHYLOGENETIC_FOREST__

#include <map>
#include <vector>

#include <Rcpp.h>

#include <phylogenetic_forest.hpp>
#include <mutation_engine.hpp>

#include "forest.hpp"

class MutationEngine;

class PhylogeneticForest : public RACES::Mutations::PhylogeneticForest
{
  std::string germline_subject;
  std::filesystem::path reference_path;

  using TimedMutationalExposure = std::map<RACES::Time, RACES::Mutations::MutationalExposure>;

  std::map<RACES::Mutations::MutationType::Type, TimedMutationalExposure> timed_exposures;

  PhylogeneticForest(const RACES::Mutations::PhylogeneticForest& orig,
                     const std::string& germline_subject,
                     const std::filesystem::path reference_path,
                     const TimedMutationalExposure& timed_SBS_exposures,
                     const TimedMutationalExposure& timed_indel_exposures);

  PhylogeneticForest(RACES::Mutations::PhylogeneticForest&& orig,
                     const std::string& germline_subject,
                     const std::filesystem::path reference_path,
                     const TimedMutationalExposure& timed_SBS_exposures,
                     const TimedMutationalExposure& timed_indel_exposures);
public:
  PhylogeneticForest();

  inline Rcpp::List get_nodes() const
  {
    return ForestCore::get_nodes(static_cast<const RACES::Mutations::PhylogeneticForest&>(*this));
  }

  const std::list<RACES::Mutants::CellId>&
  get_cell_ids_in(const std::string& sample_name) const;

  Rcpp::List get_samples_info() const;

  inline Rcpp::List get_species_info() const
  {
    return ForestCore::get_species_info(static_cast<const RACES::Mutations::PhylogeneticForest&>(*this));
  }

  inline Rcpp::List get_coalescent_cells() const
  {
    return ForestCore::get_coalescent_cells(static_cast<const RACES::Mutations::PhylogeneticForest&>(*this));
  }

  inline Rcpp::List get_coalescent_cells(const std::list<RACES::Mutants::CellId>& cell_ids) const
  {
    return ForestCore::get_coalescent_cells(static_cast<const RACES::Mutations::PhylogeneticForest&>(*this),
                                            cell_ids);
  }

  PhylogeneticForest get_subforest_for(const std::vector<std::string>& sample_names) const;

  Rcpp::List get_absolute_chromosome_positions() const;

  Rcpp::List get_germline_SIDs() const;

  Rcpp::List get_sampled_cell_SIDs() const;

  Rcpp::List get_sampled_cell_SIDs(const RACES::Mutants::CellId& cell_ids) const;

  Rcpp::List get_sampled_cell_CNAs() const;

  Rcpp::List get_sampled_cell_CNAs(const RACES::Mutants::CellId& cell_ids) const;

  Rcpp::List get_first_occurrence(const SEXP& mutation) const;

  Rcpp::List get_timed_exposures() const;

  Rcpp::List get_bulk_allelic_fragmentation(const std::string& sample_name) const;

  inline std::string get_reference_path() const
  {
    return reference_path;
  }

  void set_reference_path(const std::string reference_path);

  inline const std::string& get_germline_subject() const
  {
    return germline_subject;
  }

  inline void save(const std::string& filename) const
  {
    save(filename, false);
  }

  void save(const std::string& filename, const bool quiet) const;

  inline static PhylogeneticForest load(const std::string& filename)
  {
    return load(filename, false);
  }

  static PhylogeneticForest load(const std::string& filename, const bool quiet);

  static Rcpp::List get_SID_dataframe(const RACES::Mutations::CellGenomeMutations& cell_mutations);

  void show() const;

  friend class MutationEngine;
};

RCPP_EXPOSED_CLASS(PhylogeneticForest)

#endif // __RRACES_PHYLOGENETIC_FOREST__

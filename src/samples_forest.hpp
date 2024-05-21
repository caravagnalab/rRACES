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

#ifndef __RRACES_SAMPLES_FOREST__
#define __RRACES_SAMPLES_FOREST__

#include <vector>

#include <Rcpp.h>

#include <phylogenetic_forest.hpp>

#include "forest.hpp"

class MutationEngine;

class SamplesForest : protected Races::Mutants::DescendantsForest
{
  SamplesForest();

  Rcpp::List get_nodes(const std::vector<Races::Mutants::CellId>& cell_ids) const;

public:
  SamplesForest(const Races::Mutants::Evolutions::Simulation& simulation);

  inline Rcpp::List get_nodes() const
  {
    return ForestCore::get_nodes(static_cast<const Races::Mutants::DescendantsForest&>(*this));
  }

  Rcpp::List get_samples_info() const;

  inline Rcpp::List get_species_info() const
  {
    return ForestCore::get_species_info(static_cast<const Races::Mutants::DescendantsForest&>(*this));
  }

  inline Rcpp::List get_coalescent_cells() const
  {
    return ForestCore::get_coalescent_cells(static_cast<const Races::Mutants::DescendantsForest&>(*this));
  }

  inline Rcpp::List get_coalescent_cells(const std::list<Races::Mutants::CellId>& cell_ids) const
  {
    return ForestCore::get_coalescent_cells(static_cast<const Races::Mutants::DescendantsForest&>(*this), 
                                            cell_ids);
  }

  inline std::list<std::list<Races::Mutants::CellId>> get_sticks() const
  {
    return Races::Mutants::DescendantsForest::get_sticks();
  }

  inline std::list<std::list<Races::Mutants::CellId>> get_sticks(const double birth_threshold) const
  {
    return Races::Mutants::DescendantsForest::get_sticks(birth_threshold);
  }

  SamplesForest get_subforest_for(const std::vector<std::string>& sample_names) const;

  void save(const std::string& filename) const;

  static SamplesForest load(const std::string& filename);

  void show() const;

  friend class MutationEngine;
};



RCPP_EXPOSED_CLASS(SamplesForest)

#endif // __RRACES_SAMPLES_FOREST__

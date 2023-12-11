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


class SamplesForest : private Races::Mutants::DescendantsForest
{
  SamplesForest();

  Rcpp::List get_nodes(const std::vector<Races::Mutants::CellId>& cell_ids) const;

public:
  SamplesForest(const Races::Mutants::Evolutions::Simulation& simulation);

  Rcpp::List get_nodes() const;

  Rcpp::List get_samples_info() const;

  Rcpp::List get_species_info() const;

  Rcpp::List get_coalescent_cells() const;

  Rcpp::List get_coalescent_cells(const std::list<Races::Mutants::CellId>& cell_ids) const;

  SamplesForest get_subforest_for(const std::vector<std::string>& sample_names) const;

  void show() const;
};

#endif // __RRACES_SAMPLES_FOREST__

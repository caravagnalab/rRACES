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

#include "phylogenetic_forest.hpp"

PhylogeneticForest::PhylogeneticForest():
  Races::Mutations::PhylogeneticForest()
{}

PhylogeneticForest::PhylogeneticForest(const Races::Mutations::PhylogeneticForest& orig):
   Races::Mutations::PhylogeneticForest(orig)
{}

PhylogeneticForest PhylogeneticForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
  PhylogeneticForest forest;

  static_cast<Races::Mutations::PhylogeneticForest&>(forest) = Races::Mutations::PhylogeneticForest::get_subforest_for(sample_names);

  return forest;
}

void PhylogeneticForest::save(const std::string& filename) const
{
  Races::Archive::Binary::Out out_archive(filename);

  Races::Mutations::PhylogeneticForest::save(out_archive);
}

PhylogeneticForest PhylogeneticForest::load(const std::string& filename)
{
  PhylogeneticForest forest;

  Races::Archive::Binary::In in_archive(filename);

  static_cast<Races::Mutations::PhylogeneticForest&>(forest) = Races::Mutations::PhylogeneticForest::load(in_archive);

  return forest;
}

void PhylogeneticForest::show() const
{
  using namespace Rcpp;

  size_t num_of_leaves{0};
  for (const auto& sample: get_samples()) {
    num_of_leaves += sample.get_cell_ids().size();
  }

  Rcout << "PhylogeneticForest" << std::endl 
        << "  # of trees: " << get_roots().size() << std::endl 
        << "  # of nodes: " << num_of_nodes() << std::endl 
        << "  # of leaves: " << num_of_leaves << std::endl 
        << "  samples: {";

  std::string sep = "";
  for (const auto& sample: get_samples()) {
    Rcout << sep << "\"" << sample.get_name() << "\"";
    sep = ", ";
  }

  Rcout << "}" << std::endl;
}

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

#include "samples_forest.hpp"


SamplesForest::SamplesForest():
  Races::Mutants::DescendantsForest()
{}

// This method produces a segmentation fault and I cannot understand why
/*
SamplesForest::SamplesForest(const Simulation& simulation):
  Races::Mutants::DescendantsForest(*(simulation.sim_ptr))
{}
*/

SamplesForest::SamplesForest(const Races::Mutants::Evolutions::Simulation& simulation):
  Races::Mutants::DescendantsForest(simulation)
{}

SamplesForest SamplesForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
  SamplesForest forest;

  static_cast< Races::Mutants::DescendantsForest&>(forest) = Races::Mutants::DescendantsForest::get_subforest_for(sample_names);

  return forest;
}

void SamplesForest::save(const std::string& filename) const
{
  Races::Archive::Binary::Out out_archive(filename);

  Races::Mutants::DescendantsForest::save(out_archive);
}

SamplesForest SamplesForest::load(const std::string& filename)
{
  SamplesForest forest;

  Races::Archive::Binary::In in_archive(filename);

  static_cast<Races::Mutants::DescendantsForest&>(forest) = Races::Mutants::DescendantsForest::load(in_archive);

  return forest;
}

void SamplesForest::show() const
{
  using namespace Rcpp;

  size_t num_of_leaves{0};
  for (const auto& sample: get_samples()) {
    num_of_leaves += sample.get_cell_ids().size();
  }

  Rcout << "SamplesForest" << std::endl 
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

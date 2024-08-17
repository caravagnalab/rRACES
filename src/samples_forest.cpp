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

#include "simulation.hpp"
#include "utility.hpp"

SamplesForest::SamplesForest():
  RACES::Mutants::DescendantsForest()
{}

SamplesForest::SamplesForest(const RACES::Mutants::Evolutions::Simulation& simulation):
  RACES::Mutants::DescendantsForest(simulation)
{}

Rcpp::List SamplesForest::get_samples_info() const
{
 return SpatialSimulation::get_samples_info(get_samples());
}

SamplesForest SamplesForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
  SamplesForest forest;

  static_cast< RACES::Mutants::DescendantsForest&>(forest) = RACES::Mutants::DescendantsForest::get_subforest_for(sample_names);

  return forest;
}

void SamplesForest::save(const std::string& filename) const
{
  RACES::Archive::Binary::Out out_archive(filename);

  RACES::Mutants::DescendantsForest::save(out_archive);
}

SamplesForest SamplesForest::load(const std::string& filename)
{
  SamplesForest forest;

  if (!std::filesystem::exists(filename)) {
    throw std::domain_error("The file \"" + filename + "\" does not exist.");
  }

  if (!std::filesystem::is_regular_file(filename)) {
    throw std::domain_error("The file \"" + filename + "\" is not a regular file.");
  }

  RACES::Archive::Binary::In in_archive(filename);

  try {
    static_cast<RACES::Mutants::DescendantsForest&>(forest) = RACES::Mutants::DescendantsForest::load(in_archive);
  } catch (RACES::Archive::WrongFileFormatDescr& ex) {
    raise_error(ex, "descendants forest");
  } catch (RACES::Archive::WrongFileFormatVersion& ex) {
    raise_error(ex, "descendants forest");
  }

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

  Rcout << "}" << std::endl << std::endl;
}

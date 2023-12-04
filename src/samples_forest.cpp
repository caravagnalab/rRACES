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
  Races::Drivers::DescendantsForest()
{}

// This method produces a segmentation fault and I cannot understand why
/*
SamplesForest::SamplesForest(const Simulation& simulation):
  Races::Drivers::DescendantsForest(*(simulation.sim_ptr))
{}
*/

SamplesForest::SamplesForest(const Races::Drivers::Simulation::Simulation& simulation):
  Races::Drivers::DescendantsForest(simulation)
{}


Rcpp::List SamplesForest::get_nodes() const
{
  std::vector<Races::Drivers::CellId> cell_ids;
  cell_ids.reserve(num_of_nodes());
  for (const auto& [cell_id, cell]: get_cells()) {
    cell_ids.push_back(cell_id);
  }

  return get_nodes(cell_ids);
}

Rcpp::List SamplesForest::get_nodes(const std::vector<Races::Drivers::CellId>& cell_ids) const
{
  using namespace Rcpp;
  using namespace Races::Drivers;

  IntegerVector ids(cell_ids.size()), ancestors(cell_ids.size());
  CharacterVector genotypes(cell_ids.size()), epi_states(cell_ids.size()),
                  sample_names(cell_ids.size());
  NumericVector birth(cell_ids.size());

  size_t i{0};
  for (const auto& cell_id: cell_ids) {
    ids[i] = cell_id;
    auto cell_node = get_node(cell_id);
    if (cell_node.is_root()) {
      ancestors[i] = NA_INTEGER;
    } else {
      ancestors[i] = cell_node.parent().get_id();
    }

    genotypes[i] = cell_node.get_genotype_name();
    epi_states[i] = GenotypeProperties::signature_to_string(cell_node.get_methylation_signature());

    if (cell_node.is_leaf()) {
      sample_names[i] = cell_node.get_sample().get_name();
    } else {
      sample_names[i] = NA_STRING;
    }
    birth[i] = static_cast<const Cell&>(cell_node).get_birth_time();

    ++i;
  }

  return DataFrame::create(_["cell_id"]=ids, _["ancestor"]=ancestors,
                           _["genotype"]=genotypes, _["epistate"]=epi_states,
                           _["sample"]=sample_names, _["birth_time"]=birth);
}


template<typename SAMPLES>
Rcpp::List get_samples_info(const SAMPLES& samples)
{
  using namespace Rcpp;

  CharacterVector sample_name(samples.size());
  NumericVector time(samples.size());
  IntegerVector ymin(samples.size()), ymax(samples.size()),
                xmin(samples.size()), xmax(samples.size()),
                non_wild(samples.size());

  size_t i{0};
  for (const auto& sample : samples) {
    sample_name[i] = sample.get_name();
    time[i] = sample.get_time();
    non_wild[i] = sample.get_cell_ids().size();

    const auto& rectangle = sample.get_region();
    xmin[i] = rectangle.lower_corner.x;
    xmax[i] = rectangle.upper_corner.x;
    ymin[i] = rectangle.lower_corner.y;
    ymax[i] = rectangle.upper_corner.y;

    ++i;
  }

  return DataFrame::create(_["name"]=sample_name, _["xmin"]=xmin,
                           _["ymin"]=ymin, _["xmax"]=xmax,
                           _["ymax"]=ymax,
                           _["tumoural cells"]=non_wild,
                           _["time"]=time);
}

Rcpp::List SamplesForest::get_samples_info() const
{
  return ::get_samples_info(get_samples());
}

Rcpp::List SamplesForest::get_coalescent_cells() const
{
  auto coalencent_ids = Races::Drivers::DescendantsForest::get_coalescent_cells();

  return get_nodes(coalencent_ids);
}

Rcpp::List SamplesForest::get_coalescent_cells(const std::list<Races::Drivers::CellId>& cell_ids) const
{
  auto coalencent_ids = Races::Drivers::DescendantsForest::get_coalescent_cells(cell_ids);

  return get_nodes(coalencent_ids);
}

SamplesForest SamplesForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
  SamplesForest forest;

  static_cast< Races::Drivers::DescendantsForest&>(forest) = Races::Drivers::DescendantsForest::get_subforest_for(sample_names);

  return forest;
}

Rcpp::List SamplesForest::get_species_info() const
{
  using namespace Rcpp;

  size_t num_of_rows = get_species_data().size();

  CharacterVector genotype_names(num_of_rows), epi_states(num_of_rows);

  using namespace Races::Drivers;

  size_t i{0};
  for (const auto& [species_id, species_data]: get_species_data()) {
    genotype_names[i] = get_genotype_name(species_data.genotype_id);
    epi_states[i] = GenotypeProperties::signature_to_string(species_data.signature);

    ++i;
  }

  return DataFrame::create(_["genotype"]=genotype_names, _["epistate"]=epi_states);
}

void SamplesForest::show() const
{
  using namespace Rcpp;

  size_t num_of_leaves{0};
  for (const auto& sample: get_samples()) {
    num_of_leaves += sample.get_cell_ids().size();
  }

  Rcout << "SamplesForest(# of trees: " << get_roots().size()
        << ", # of nodes: " << num_of_nodes()
        << ", # of leaves: " << num_of_leaves
        << ", samples: {";

  std::string sep = "";
  for (const auto& sample: get_samples()) {
    Rcout << sep << "\"" << sample.get_name() << "\"";
    sep = ", ";
  }

  Rcout << "})" << std::endl;
}

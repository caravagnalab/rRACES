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

#ifndef __RRACES_FOREST__
#define __RRACES_FOREST__

#include <vector>

#include <Rcpp.h>

#include <mutant_properties.hpp>

struct ForestCore
{
  template<typename CPP_FOREST>
  static Rcpp::List get_nodes(const CPP_FOREST& forest,
                              const std::vector<RACES::Mutants::CellId>& cell_ids)
  {
    using namespace Rcpp;
    using namespace RACES::Mutants;

    IntegerVector ids(cell_ids.size()), ancestors(cell_ids.size());
    CharacterVector mutants(cell_ids.size()), epi_states(cell_ids.size()),
                    sample_names(cell_ids.size());
    NumericVector birth(cell_ids.size());

    size_t i{0};
    for (const auto& cell_id: cell_ids) {
      ids[i] = cell_id;
      auto cell_node = forest.get_node(cell_id);
      if (cell_node.is_root()) {
        ancestors[i] = NA_INTEGER;
      } else {
        ancestors[i] = cell_node.parent().get_id();
      }

      mutants[i] = cell_node.get_mutant_name();
      epi_states[i] = MutantProperties::signature_to_string(cell_node.get_methylation_signature());

      if (cell_node.is_leaf()) {
        sample_names[i] = cell_node.get_sample().get_name();
      } else {
        sample_names[i] = NA_STRING;
      }
      birth[i] = static_cast<const Cell&>(cell_node).get_birth_time();

      ++i;
    }

    return DataFrame::create(_["cell_id"]=ids, _["ancestor"]=ancestors,
                             _["mutant"]=mutants, _["epistate"]=epi_states,
                             _["sample"]=sample_names, _["birth_time"]=birth);
  }

  template<typename CPP_FOREST>
  static Rcpp::List get_nodes(const CPP_FOREST& forest)
  {
    std::vector<RACES::Mutants::CellId> cell_ids;
    cell_ids.reserve(forest.num_of_nodes());
    for (const auto& [cell_id, cell]: forest.get_cells()) {
      cell_ids.push_back(cell_id);
    }

    return ForestCore::get_nodes<CPP_FOREST>(forest, cell_ids);
  }

  template<typename CPP_FOREST>
  static Rcpp::List get_species_info(const CPP_FOREST& forest)
  {
    using namespace Rcpp;

    size_t num_of_rows = forest.get_species_data().size();

    CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);

    using namespace RACES::Mutants;

    size_t i{0};
    for (const auto& [species_id, species_data]: forest.get_species_data()) {
      mutant_names[i] = forest.get_mutant_name(species_data.mutant_id);
      epi_states[i] = MutantProperties::signature_to_string(species_data.signature);

      ++i;
    }

    return DataFrame::create(_["mutant"]=mutant_names, _["epistate"]=epi_states);
  }

  template<typename CPP_FOREST>
  static Rcpp::List get_coalescent_cells(const CPP_FOREST& forest)
  {
    auto coalencent_ids = forest.get_coalescent_cells();

    return ForestCore::get_nodes<CPP_FOREST>(forest, coalencent_ids);
  }

  template<typename CPP_FOREST>
  static Rcpp::List get_coalescent_cells(const CPP_FOREST& forest,
                                        const std::list<RACES::Mutants::CellId>& cell_ids)
  {
    auto coalencent_ids = forest.get_coalescent_cells(cell_ids);

    return ForestCore::get_nodes<CPP_FOREST>(forest, coalencent_ids);
  }
};

#endif // __RRACES_FOREST__

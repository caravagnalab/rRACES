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

#include "sampled_cell.hpp"

#include "phylogenetic_forest.hpp"

SampledCell::SampledCell(const Races::Mutations::PhylogeneticForest& forest,
           const Races::Mutants::CellId& cell_id):
           Races::Mutations::PhylogeneticForest::const_node(&forest, cell_id)
{}

Rcpp::List SampledCell::mutations() const
{
    return PhylogeneticForest::get_SID_dataframe(cell_mutations());
}
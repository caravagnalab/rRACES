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

#ifndef __RRACES_SAMPLED_CELL__
#define __RRACES_SAMPLED_CELL__

#include <string>

#include <phylogenetic_forest.hpp>

#include <Rcpp.h>

class SampledCell : private Races::Mutations::PhylogeneticForest::const_node
{
public:
    SampledCell(const Races::Mutations::PhylogeneticForest& forest,
         const Races::Mutants::CellId& cell_id);

    inline std::string epistate() const
    {
        using namespace Races::Mutants;

        return MutantProperties::signature_to_string(get_methylation_signature());
    }

    inline std::string mutant() const
    {
        return get_mutant_name();
    }

    inline std::string species() const
    {
        return mutant() + epistate();
    }

    inline const Races::Time& birth_time() const
    {
        return get_birth_time();
    }

    Rcpp::List mutations() const;
};

RCPP_EXPOSED_CLASS(SampledCell)

#endif // __RRACES_SAMPLED_CELL__

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

#ifndef __RRACES_WG_DOUBLING__
#define __RRACES_WG_DOUBLING__

#include <string>

#include <cna.hpp>

#include <Rcpp.h>


class WholeGenomeDoubling
{
public:
    inline WholeGenomeDoubling()
    {}

    inline void show() const
    {
        using namespace Rcpp;

        Rcout << "Whole Genome Doubling";
    }
};

RCPP_EXPOSED_CLASS(WholeGenomeDoubling)

#endif // __RRACES_WG_DOUBLING__

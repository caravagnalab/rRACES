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

#ifndef __RRACES_TISSUE_RECTANGLE__
#define __RRACES_TISSUE_RECTANGLE__

#include <cstdint>
#include <vector>

#include <Rcpp.h>

#include <position_set.hpp>
#include <position.hpp>


class TissueRectangle : public RACES::Mutants::RectangleSet
{
public:
  TissueRectangle();

  TissueRectangle(const RACES::Mutants::Evolutions::PositionInTissue& lower_corner, 
                  const RACES::Mutants::Evolutions::PositionInTissue& upper_corner);

  TissueRectangle(const RACES::Mutants::Evolutions::PositionInTissue& lower_corner, 
                  const RACES::Mutants::Evolutions::AxisSize& x_size, 
                  const RACES::Mutants::Evolutions::AxisSize& y_size);

  TissueRectangle(const std::vector<uint16_t>& lower_corner, 
                  const std::vector<uint16_t>& upper_corner);

  TissueRectangle(const std::vector<uint16_t>& lower_corner, 
                  const RACES::Mutants::Evolutions::AxisSize& x_size, 
                  const RACES::Mutants::Evolutions::AxisSize& y_size);

  Rcpp::IntegerVector get_lower_corner() const;

  Rcpp::IntegerVector get_upper_corner() const;

  void show() const;
};

RCPP_EXPOSED_CLASS(TissueRectangle)

#endif // __RRACES_TISSUE_RECTANGLE__

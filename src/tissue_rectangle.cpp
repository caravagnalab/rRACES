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

#include "tissue_rectangle.hpp"

Rcpp::IntegerVector TissueRectangle::get_lower_corner() const
{
  Rcpp::IntegerVector lcorner(2);

  lcorner[0] = lower_corner.x;
  lcorner[1] = lower_corner.y;

  return lcorner;
}

Rcpp::IntegerVector TissueRectangle::get_upper_corner() const
{
  Rcpp::IntegerVector ucorner(2);

  ucorner[0] = upper_corner.x;
  ucorner[1] = upper_corner.y;

  return ucorner;
}

void TissueRectangle::show() const
{
  using namespace Rcpp;

  Rcout << "TissueRectangle("
        << "(" << lower_corner.x <<"," << lower_corner.y << "),"
        << "(" << upper_corner.x <<"," << upper_corner.y << "))" << std::endl;
}

TissueRectangle::TissueRectangle():
  RectangleSet()
{}

TissueRectangle::TissueRectangle(const RACES::Mutants::Evolutions::PositionInTissue& lower_corner, 
                                 const RACES::Mutants::Evolutions::PositionInTissue& upper_corner):
  RectangleSet(lower_corner, upper_corner)
{}

TissueRectangle::TissueRectangle(const RACES::Mutants::Evolutions::PositionInTissue& lower_corner, 
                                 const RACES::Mutants::Evolutions::AxisSize& x_size, 
                                 const RACES::Mutants::Evolutions::AxisSize& y_size):
  RectangleSet(lower_corner, x_size, y_size)
{}

TissueRectangle::TissueRectangle(const std::vector<uint16_t>& lower_corner, const std::vector<uint16_t>& upper_corner):
  TissueRectangle(RACES::Mutants::Evolutions::PositionInTissue{lower_corner[0],lower_corner[1]},
                  RACES::Mutants::Evolutions::PositionInTissue{upper_corner[0], upper_corner[1]})
{}

TissueRectangle::TissueRectangle(const std::vector<uint16_t>& lower_corner, 
                                 const RACES::Mutants::Evolutions::AxisSize& x_size, 
                                 const RACES::Mutants::Evolutions::AxisSize& y_size):
  TissueRectangle(RACES::Mutants::Evolutions::PositionInTissue{lower_corner[0],lower_corner[1]},
                  x_size, y_size)
{}


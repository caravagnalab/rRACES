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

//' @name TissueRectangle
//' @title A rectangle in the tissue
//' @field get_lower_corner Get the rectangle lower corner
//' @field get_upper_corner Get the rectangle upper corner
class TissueRectangle : public Races::Drivers::RectangleSet
{
public:
  TissueRectangle(const Races::Drivers::Simulation::PositionInTissue& lower_corner, 
                  const Races::Drivers::Simulation::PositionInTissue& upper_corner);

  TissueRectangle(const Races::Drivers::Simulation::PositionInTissue& lower_corner, 
                  const Races::Drivers::Simulation::AxisSize& x_size, 
                  const Races::Drivers::Simulation::AxisSize& y_size);

  TissueRectangle(const std::vector<uint16_t>& lower_corner, 
                  const std::vector<uint16_t>& upper_corner);

  TissueRectangle(const std::vector<uint16_t>& lower_corner, 
                  const Races::Drivers::Simulation::AxisSize& x_size, 
                  const Races::Drivers::Simulation::AxisSize& y_size);

  Rcpp::IntegerVector get_lower_corner() const;

  Rcpp::IntegerVector get_upper_corner() const;

  void show() const;
};

//' @name TissueRectangle$new
//' @title Build a new rectangle of tissue.
//' @examples
//' # build the rectangle [500,550]x[450,475]
//' rect <- new(TissueRectangle, c(500, 450), c(550, 475))
//'
//' rect
//'
//' # build the rectangle [500,550]x[450,475]
//' rect <- new(TissueRectangle, c(500, 450), 50, 25)
//'
//' rect

//' @name TissueRectangle$lower_corner
//' @title The lower corner of the tissue rectangle.
//' @examples
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the simulation death activation level
//' rect$lower_corner

//' @name TissueRectangle$upper_corner
//' @title The lower corner of the tissue rectangle.
//' @examples
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the simulation death activation level
//' rect$upper_corner

#endif // __RRACES_TISSUE_RECTANGLE__

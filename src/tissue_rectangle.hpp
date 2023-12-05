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

#endif // __RRACES_TISSUE_RECTANGLE__

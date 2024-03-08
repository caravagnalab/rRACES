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

#ifndef __RRACES_SEQ_SIMULATION__
#define __RRACES_SEQ_SIMULATION__

#include <Rcpp.h>

#include <read_simulator.hpp>

#include "phylogenetic_forest.hpp"

Rcpp::List  simulate_seq(const PhylogeneticForest& forest, const double& coverage, 
                         const int& read_size, const int& insert_size,
                         const std::string& output_dir, const bool& write_SAM,
                         const bool& FACS, const double& purity, const int& rnd_seed);

#endif // __RRACES_SEQ_SIMULATION__

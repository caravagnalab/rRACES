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

#ifndef __RRACES_SETUP_MUTATION_ENGINE__
#define __RRACES_SETUP_MUTATION_ENGINE__

#include <string>

#include <Rcpp.h>


Rcpp::List get_mutation_engine_supported_codes();

void setup_mutation_engine_by_code(const std::string& setup_code,
                                     const size_t& context_sampling_rate);

void setup_mutation_engine_by_code(const std::string& setup_code);

void setup_mutation_engine(const std::string& directory,
                           const std::string& reference_url,
                           const std::string& SBS_url,
                           const size_t& context_sampling);

void setup_mutation_engine(const std::string& directory,
                           const std::string& reference_url,
                           const std::string& SBS_url);

#endif // __RRACES_SETUP_MUTATION_ENGINE__

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

#include "sequencers.hpp"

ErrorlessIlluminaSequencer::ErrorlessIlluminaSequencer():
    Races::Sequencers::Illumina::ErrorLessSequencer()
{}

void ErrorlessIlluminaSequencer::show() const
{
    Rcpp::Rcout << get_model_name() << " (platform: \""
                << get_platform_name() << "\")" << std::endl;
}

BasicIlluminaSequencer::BasicIlluminaSequencer(const double error_rate, const int seed):
    Races::Sequencers::Illumina::BasicSequencer<>(error_rate,seed)
{}

void BasicIlluminaSequencer::show() const
{
    Rcpp::Rcout << get_model_name() << " (platform: \""
                << get_platform_name() << "\" error rate: " 
                << std::to_string(get_error_rate()) << ")" << std::endl;
}


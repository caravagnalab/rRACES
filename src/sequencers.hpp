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

#ifndef __RRACES_SEQUENCERS__
#define __RRACES_SEQUENCERS__

#include <Rcpp.h>

#include <sequencer.hpp>

class ErrorlessIlluminaSequencer : public RACES::Sequencers::Illumina::ErrorLessSequencer
{
public:
    ErrorlessIlluminaSequencer();

    void show() const;
};

class BasicIlluminaSequencer : public RACES::Sequencers::Illumina::BasicSequencer<>
{
public:
    BasicIlluminaSequencer(const double error_rate, const SEXP seed=R_NilValue);

    void show() const;
};

RCPP_EXPOSED_CLASS(ErrorlessIlluminaSequencer)
RCPP_EXPOSED_CLASS(BasicIlluminaSequencer)

#endif // __RRACES_SEQUENCERS__

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

    static ErrorlessIlluminaSequencer build_sequencer();
};

class BasicIlluminaSequencer
{
    std::shared_ptr<RACES::Sequencers::Illumina::BasicSequencer<RACES::Sequencers::Illumina::BasicQualityScoreModel>> random_score_seq;
    std::shared_ptr<RACES::Sequencers::Illumina::BasicSequencer<RACES::Sequencers::ConstantQualityScoreModel>> constant_score_seq;
public:

    BasicIlluminaSequencer(const double error_rate, const bool& random_quality_scores=true,
                           const int seed=0);

    void show() const;

    const double& get_error_rate() const;

    inline bool producing_random_scores() const
    {
        return random_score_seq != nullptr;
    }

    template<typename QUALITY_SCORE_MODEL>
    RACES::Sequencers::Illumina::BasicSequencer<QUALITY_SCORE_MODEL>& basic_sequencer() const
    {
        if constexpr (std::is_base_of_v<QUALITY_SCORE_MODEL, RACES::Sequencers::Illumina::BasicQualityScoreModel>) {
            return *random_score_seq;
        }
        if constexpr (std::is_base_of_v<QUALITY_SCORE_MODEL, RACES::Sequencers::ConstantQualityScoreModel>) {
            return *constant_score_seq;
        }

        throw std::domain_error("Unsupported quality score model.");
    }

    static BasicIlluminaSequencer build_sequencer(const SEXP error_rate,
                                                  const SEXP random_quality_scores,
                                                  const SEXP seed);
};

RCPP_EXPOSED_CLASS(ErrorlessIlluminaSequencer)
RCPP_EXPOSED_CLASS(BasicIlluminaSequencer)

#endif // __RRACES_SEQUENCERS__

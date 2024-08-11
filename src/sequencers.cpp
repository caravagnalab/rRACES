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

#include "utility.hpp"

ErrorlessIlluminaSequencer::ErrorlessIlluminaSequencer():
    RACES::Sequencers::Illumina::ErrorLessSequencer()
{}

void ErrorlessIlluminaSequencer::show() const
{
    Rcpp::Rcout << get_model_name() << " (platform: \""
                << get_platform_name() << "\")" << std::endl;
}

ErrorlessIlluminaSequencer
ErrorlessIlluminaSequencer::build_sequencer()
{
    return ErrorlessIlluminaSequencer();
}

BasicIlluminaSequencer::BasicIlluminaSequencer(const double error_rate,
                                               const bool& random_quality_scores,
                                               const int seed):
    random_score_seq(nullptr), constant_score_seq(nullptr)
{
    using namespace RACES::Sequencers;
    if (random_quality_scores) {
        random_score_seq = std::make_shared<Illumina::BasicSequencer<>>(error_rate, seed);
    } else {
        constant_score_seq = std::make_shared<Illumina::BasicSequencer<ConstantQualityScoreModel>>(error_rate, seed);
    }
}

template<typename QUALITY_SCORE_MODEL>
void show_sequencer(const RACES::Sequencers::Illumina::BasicSequencer<QUALITY_SCORE_MODEL>& seq)
{
    Rcpp::Rcout << seq.get_model_name() << " (platform: \""
                << seq.get_platform_name() << "\" error rate: " 
                << std::to_string(seq.get_error_rate());

    using namespace RACES::Sequencers;

    if constexpr (std::is_base_of_v<QUALITY_SCORE_MODEL,Illumina::BasicQualityScoreModel> ) {            
        Rcpp::Rcout << " random quality scores";
    }

    if constexpr (std::is_base_of_v<QUALITY_SCORE_MODEL,ConstantQualityScoreModel> ) {            
        Rcpp::Rcout << " constant quality scores" << std::endl;
    }

    Rcpp::Rcout << ")" << std::endl;
}

void BasicIlluminaSequencer::show() const
{
    if (producing_random_scores()) {
        show_sequencer(*random_score_seq);
    } else {
        show_sequencer(*constant_score_seq);
    }
}

const double& BasicIlluminaSequencer::get_error_rate() const
{
    if (producing_random_scores()) {
        return random_score_seq->get_error_rate();
    } else {
        return constant_score_seq->get_error_rate();
    }     
}

BasicIlluminaSequencer
BasicIlluminaSequencer::build_sequencer(const SEXP error_rate,
                                        const SEXP random_quality_scores,
                                        const SEXP seed)
{
    using namespace Rcpp;

    if (TYPEOF(error_rate) != INTSXP && TYPEOF(error_rate) != REALSXP) {
        throw std::domain_error("The parameter \"error_rate\""
                                " must be a positive real number.");
    }

    const auto c_error_rate = as<double>(error_rate); 

    if (c_error_rate<0) {
        throw std::domain_error("The parameter \"error_rate\""
                                " must be a positive real number.");
    }

    if (TYPEOF(random_quality_scores) != LGLSXP) {
        throw std::domain_error("The parameter \"random_quality_scores\""
                                " must be a Boolean value.");
    }

    const auto c_random_quality_scores = as<bool>(random_quality_scores); 

    auto c_seed = get_random_seed<int>(seed);

    return {c_error_rate, c_random_quality_scores, c_seed};
}
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

#include "snv.hpp"

SNV::SNV(const Races::Mutations::ChromosomeId& chromosome_id,
         const Races::Mutations::ChrPosition& chromosomic_position,
         const char& ref_base, const char& alt_base, const std::string& cause):
    Races::Mutations::SNV(chromosome_id, chromosomic_position, ref_base,
                          alt_base, cause)
{}

SNV::SNV()
{}

SEXP SNV::get_cause() const
{
    Rcpp::StringVector cause_v(1);

    if (cause=="") {
        cause_v[0] = NA_STRING;
    } else {
        cause_v[0] = cause;
    }

    return cause_v;
}

Rcpp::List SNV::get_dataframe() const
{
    using namespace Rcpp;
    using namespace Races::Mutations;

    return DataFrame::create(_["chromosome"]=get_chromosome(),
                             _["pos_in_chr"]=position,
                             _["ref"]=get_ref_base(),
                             _["alt"]=get_alt_base(),
                             _["cause"]=get_cause());
}

void SNV::show() const
{
    using namespace Rcpp;

    Rcout << "SNV(chromosome: "<< get_chromosome()
          << ", pos_in_chr: " << static_cast<size_t>(position)
          << ", ref: " << ref_base
          << ", alt: " << alt_base;

    if (cause!="") {
        Rcout << ", cause: \"" << cause << "\"";
    }
    Rcout << ")" << std::endl;
}

SNV SNV::build_SNV(const SEXP chromosome_name,
                   const SEXP position_in_chromosome, 
                   const SEXP alt_base, const SEXP ref_base, 
                   const SEXP cause)
{
    using namespace Rcpp;
    using namespace Races::Mutations;

    auto chr_name = as<std::string>(chromosome_name);
    auto chr_id = GenomicPosition::stochr(chr_name);

    auto pos = Rcpp::as<long int>(position_in_chromosome);
    if (pos < 0) {
        throw std::domain_error("Position in chromosome must be a "
                                "non-negative number");
    }

    auto ref_base_str = Rcpp::as<std::string>(ref_base);
    if (ref_base_str.size() != 1) {
        throw std::domain_error("The reference base must be a single "
                                "nucleotide.");
    }

    auto alt_base_str = Rcpp::as<std::string>(alt_base);
    if (alt_base_str.size() != 1) {
        throw std::domain_error("The altered base must be a single "
                                "nucleotide.");
    }

    auto cause_str = Rcpp::as<std::string>(cause);
    return SNV(chr_id, static_cast<Races::Mutations::ChrPosition>(pos),
               ref_base_str[0], alt_base_str[0], cause_str);
}
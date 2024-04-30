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

#include "sid.hpp"

#include "utility.hpp"

SID::SID(const Races::Mutations::ChromosomeId& chromosome_id,
         const Races::Mutations::ChrPosition& chromosomic_position,
         const Races::Mutations::AlleleId allele_id,
         const char& ref, const char& alt, const std::string& cause):
    Races::Mutations::MutationSpec<Races::Mutations::SID>(allele_id, chromosome_id,
                                                          chromosomic_position,
                                                          ref, alt, cause)
{}

std::string alleletostr(const Races::Mutations::AlleleId& allele_id)
{
    if (allele_id == RANDOM_ALLELE) {
        return "random";
    }

    return std::to_string(allele_id);
}

SEXP wrap_allele(const Races::Mutations::AlleleId& allele_id)
{
    Rcpp::StringVector allele_v(1);

    allele_v[0] = alleletostr(allele_id);

    return allele_v;
}

SID::SID()
{}

SEXP SID::get_cause() const
{
    Rcpp::StringVector cause_v(1);

    if (cause=="") {
        cause_v[0] = NA_STRING;
    } else {
        cause_v[0] = cause;
    }

    return cause_v;
}

Rcpp::List SID::get_dataframe() const
{
    using namespace Rcpp;
    using namespace Races::Mutations;

    return DataFrame::create(_["chr"]=get_chromosome(),
                             _["chr_pos"]=position,
                             _["allele"]=wrap_allele(allele_id),
                             _["ref"]=get_ref(),
                             _["alt"]=get_alt(),
                             _["type"]=(is_SNV()?"SNV":"indel"),
                             _["cause"]=get_cause());
}

void SID::show() const
{
    using namespace Rcpp;

    if (is_SNV()) {
        Rcout << "SNV";
    } else {
        Rcout << "indel";
    }

    Rcout << "(chr: "<< get_chromosome()
          << ", chr_pos: " << static_cast<size_t>(position)
          << ", allele: " << alleletostr(allele_id)
          << ", ref: " << (ref.size()==0?"-":ref)
          << ", alt: " << (alt.size()==0?"-":alt);

    if (cause!="") {
        Rcout << ", cause: \"" << cause << "\"";
    }
    Rcout << ")" << std::endl;
}

SID SID::build_SNV(const SEXP chromosome_name,
                   const SEXP position_in_chromosome,
                   const SEXP alt_base, const SEXP ref_base,
                   const SEXP allele_id, const SEXP cause)
{
    SID mutation = build_SID(chromosome_name, position_in_chromosome,
                             ref_base, alt_base, allele_id, cause);

    if (mutation.ref.size() != 1) {
        throw std::domain_error("The reference base must be a single "
                                "nucleotide.");
    }

    if (mutation.alt.size() != 1) {
        throw std::domain_error("The altered base must be a single "
                                "nucleotide.");
    }

    return mutation;
}

SID SID::build_SID(const SEXP chromosome_name,
                   const SEXP position_in_chromosome,
                   const SEXP ref_base, const SEXP alt_base,
                   const SEXP allele_id, const SEXP cause)
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
    auto alt_base_str = Rcpp::as<std::string>(alt_base);

    auto cause_str = Rcpp::as<std::string>(cause);
    return SID(chr_id, static_cast<Races::Mutations::ChrPosition>(pos),
               get_allele_id(allele_id, "allele"), ref_base_str[0],
               alt_base_str[0], cause_str);
}
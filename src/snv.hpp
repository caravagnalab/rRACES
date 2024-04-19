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

#ifndef __RRACES_SNV__
#define __RRACES_SNV__

#include <string>

#include <snv.hpp>
#include <mutation_spec.hpp>

#include <Rcpp.h>

class SNV : public Races::Mutations::MutationSpec<Races::Mutations::SNV>
{

    SNV(const Races::Mutations::ChromosomeId& chromosome_id,
        const Races::Mutations::ChrPosition& chromosomic_position,
        const Races::Mutations::AlleleId allele_id,
        const char& ref_base, const char& alt_base,
        const std::string& cause="");

public:
    SNV();

    inline std::string get_chromosome() const
    {
        return Races::Mutations::GenomicPosition::chrtos(chr_id);
    }

    inline const Races::Mutations::ChrPosition& get_position_in_chromosome() const
    {
        return position;
    }

    inline std::string get_ref_base() const
    {
        return std::string(1,ref_base);
    }

    inline std::string get_alt_base() const
    {
        return std::string(1,alt_base);
    }

    SEXP get_cause() const;

    Rcpp::List get_dataframe() const;

    void show() const;

    static
    SNV build_SNV(const SEXP chromosome_name,
                  const SEXP position_in_chromosome,
                  const SEXP alt_base, const SEXP ref_base,
                  const SEXP allele_id, const SEXP cause);
};

RCPP_EXPOSED_CLASS(SNV);

#endif // __RRACES_SNV__

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

#ifndef __RRACES_SNV__
#define __RRACES_SNV__

#include <string>

#include <snv.hpp>

#include <Rcpp.h>


class SNV : public Races::Mutations::SNV
{

    SNV(const Races::Mutations::ChromosomeId& chromosome_id,
        const Races::Mutations::ChrPosition& chromosomic_position,
        const char* context, const char& mutated_base,
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

    inline std::string get_context() const
    {
        return context.get_sequence();
    }

    inline std::string get_mutated_base() const
    {
        return std::string(1,mutated_base);
    }

    inline SEXP get_cause() const
    {
        if (cause!="") {
            Rcpp::CharacterVector cause_v(1);

            cause_v[0]=cause;

            return cause_v;
        }

        return NA_STRING;
    }

    Rcpp::List get_dataframe() const;

    void show() const;

    static
    SNV build_SNV(const SEXP chromosome_name,
                  const SEXP position_in_chromosome,
                  const SEXP context, const SEXP mutated_base,
                  const SEXP cause);
};

RCPP_EXPOSED_CLASS(SNV);

#endif // __RRACES_SNV__

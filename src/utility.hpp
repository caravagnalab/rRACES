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

#ifndef __RRACES_UTILITY__
#define __RRACES_UTILITY__

#include <string>
#include <filesystem>

#include <Rcpp.h>

#include <allele.hpp>

std::filesystem::path get_tmp_dir_path(const std::string& base_name="rRACES");

Races::Mutations::AlleleId get_allele_id(const SEXP allele_id,
                                         const std::string& parameter_name);

std::string ordinal_suffix(const size_t& ord);

inline std::string ordtostr(const size_t ord)
{
    return std::to_string(ord) + ordinal_suffix(ord);
}

#endif // __RRACES_UTILITY__

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

#include <iostream>

#include <sstream>
#include <random>
#include <cstdint>

#include <utils.hpp>


#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#include <unistd.h>
#include <pwd.h>
#elif defined (__WIN32__) || defined(__WIN64__)
#include <Windows.h>
#endif

#include "utility.hpp"

std::string get_user_name()
{
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))

    auto userid = getuid();
    auto pwd = getpwuid(userid);
    return pwd->pw_name;

#elif defined (__WIN32__) || defined(__WIN64__)

    const int MAX_LEN = 100;
    char szBuffer[MAX_LEN];
    DWORD len = MAX_LEN;
    if (GetUserName(szBuffer, &len)) {
        return szBuffer;
    }

    return "unknown";

#else

    return "unknown";

#endif
}

template<typename T, typename std::enable_if_t<std::is_integral_v<T>, bool> = true>
std::string int2hex(const T value)
{
  std::ostringstream oss;

  oss << std::setfill('0') << std::setw(8) << std::hex << value;

  return oss.str();
}

std::filesystem::path get_tmp_dir_path(const std::string& base_name)
{
  int i{0};

  std::mt19937 r_gen((unsigned) time(NULL));

  std::uniform_int_distribution<uint32_t> dist;

  auto base_path = to_string(std::filesystem::temp_directory_path()/
                               (base_name + "_" + get_user_name()));

  std::filesystem::path output_path = base_path + "_" + int2hex(dist(r_gen));

  while (std::filesystem::exists(output_path)) {
    output_path = base_path + "_" + int2hex(dist(r_gen)+(++i));
  }

  return output_path;
}

Races::Mutations::AlleleId
get_allele_id(const SEXP allele_id, const std::string& parameter_name)
{
    switch (TYPEOF(allele_id)) {
        case REALSXP:
        {
            long int allele_l;
            try {
                allele_l = Rcpp::as<long int>(allele_id);
            } catch (std::invalid_argument& ex) {
                allele_l = -1;
            }

            if (allele_l >= 0) {
                return allele_l;
            }
            break;
        }
        case NILSXP:
        {
            return RANDOM_ALLELE;
        }
        default:
            break;
    }

    throw std::domain_error("The parameter \"" + parameter_name
                            + "\" must be either a "
                            + "non-negative number or NILL.");
}
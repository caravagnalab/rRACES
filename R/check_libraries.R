## This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
## Copyright (C) 2023 - Alberto Casagrande <alberto.casagrande@uniud.it>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

check_library <- function(library_name) {
  return(nzchar(system.file(package = library_name)))
}

check_libraries <- function(library_names) {
  index <- 1

  while (length(library_names) >= index) {
    if (!check_library(library_names[index])) {
      return(FALSE)
    }

    index <- index + 1
  }

  return(TRUE)
}

required_packages <- function(function_name, dependencies) {

    for (dep in dependencies) {
        if (!requireNamespace(dep, quietly = TRUE)) {
            stop(paste0("\"", function_name, "\" requires \"", dep,
                        "\". Please install it to call \"", function_name,"\"."));
        }
    }
}

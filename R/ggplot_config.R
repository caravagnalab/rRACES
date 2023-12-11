## This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
## Copyright (C) 2023 - Giulio Caravagna <gcaravagna@units.it>
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

my_theme <- function() {
  ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

get_colors_for <- function(df, column_name, pal_name) {
  # Unpaired
  colors <- NULL
  if (nrow(df) > 0) {
    num_of_values <- df[, column_name] %>% length()

    if (num_of_values < 3) {
      colors <- RColorBrewer::brewer.pal(3, pal_name)

      colors <- colors[seq_len(num_of_values)]
    } else {
      colors <- RColorBrewer::brewer.pal(num_of_values, pal_name)
    }

    names(colors) <- df[, column_name]
  }

  return(colors)
}

get_species_colors <- function(species) {

  paired_species <- species %>%
    dplyr::filter(!is.na(.data$epistate)) %>%
    dplyr::mutate(
      species = paste0(.data$mutant, .data$epistate)
    ) %>%
    dplyr::arrange(.data$mutant)

  unpaired_species <- species %>%
    dplyr::filter(is.na(.data$epistate)) %>%
    dplyr::mutate(
      species = paste0(.data$mutant, .data$epistate)
    ) %>%
    dplyr::arrange(.data$mutant)

  # Unpaired
  unp_colour <- get_colors_for(unpaired_species, "mutant", "Dark2")

  # Paired
  pai_colour <- get_colors_for(paired_species, "species", "Paired")

  return(c(pai_colour, unp_colour))
}





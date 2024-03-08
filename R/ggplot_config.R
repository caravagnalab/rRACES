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

#' @import Polychrome
get_signatures_colors <-  function(seed = 1999) {
  set.seed(seed)
  sigs = c('SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS10a', 'SBS10b', 'SBS11', 'SBS13', 'SBS14', 'SBS15', 'SBS17b', 'SBS18', 'SBS20', 'SBS21', 'SBS22a', 'SBS22b', 'SBS24', 'SBS25', 'SBS26', 'SBS28', 'SBS30', 'SBS31', 'SBS35', 'SBS36', 'SBS40a', 'SBS40b', 'SBS40c', 'SBS42', 'SBS87', 'SBS88', 'SBS90', 'SBS96', 'SBS97', 'SBS98', 'SBS99', 'SBS5', 'DBS1', 'DBS2', 'DBS3', 'DBS4', 'DBS5', 'DBS6', 'DBS7', 'DBS8', 'DBS9', 'DBS10', 'DBS11', 'DBS12', 'DBS13', 'DBS14', 'DBS15', 'DBS16', 'DBS17', 'DBS18', 'DBS19', 'DBS20', 'CN1', 'CN2', 'CN3', 'CN4', 'CN5', 'CN6', 'CN7', 'CN8', 'CN9', 'CN10', 'CN11', 'CN12', 'CN13', 'CN14', 'CN15', 'CN16', 'CN17', 'CN18', 'CN19', 'CN20', 'CN21', 'CN22', 'CN23', 'CN24', 'ID1', 'ID2', 'ID3', 'ID4', 'ID5', 'ID6', 'ID7', 'ID8', 'ID9', 'ID10', 'ID11', 'ID12', 'ID13', 'ID14', 'ID15', 'ID16', 'ID17', 'ID18', 'SBS8', 'SBS9', 'SBS10c', 'SBS10d', 'SBS12', 'SBS16', 'SBS17a', 'SBS19', 'SBS23', 'SBS27', 'SBS29', 'SBS32', 'SBS33', 'SBS34', 'SBS37', 'SBS38', 'SBS39', 'SBS41', 'SBS43', 'SBS44', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'SBS84', 'SBS85', 'SBS86', 'SBS89', 'SBS91', 'SBS92', 'SBS93', 'SBS94', 'SBS95') 
  return(Polychrome::createPalette(length(sigs), c("#6B8E23","#4169E1"), 
                                   M = 1000,
                                   target = "normal", range = c(15,80)) %>%
           setNames(sigs))
}




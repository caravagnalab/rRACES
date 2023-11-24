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

#' Plot a tissue
#'
#' @description
#' Plots cells distribution over a tissue highlighting species by color.
#' To facilitate the plot and avoid excessive number of cells, for instance,
#' when a simulation deals with millions of cells, the plot draws a
#' hexagonal heatmap of 2D bins.
#'
#' @param simulation A simulation object.
#' @param num_of_bins The number of bins (optional, default: 100).
#'
#' @return An editable ggplot plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation)
#' sim$add_genotype(name = "A",
#'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
#'                  growth_rates = c("+" = 0.2, "-" = 0.08),
#'                  death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#' plot_tissue(sim)
#'
#' # delete the simulation dump directory
#' unlink(sim$get_name(), recursive=TRUE)
plot_tissue <- function(simulation, num_of_bins = 100) {
  stopifnot(inherits(simulation, "Rcpp_Simulation"))

  # Get the cells in the tissue at current simulation time
  cells <- simulation$get_cells() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(species = paste0(.data$genotype, .data$epistate))

  counts <- simulation$get_counts()

  time <- simulation$get_clock() %>% round(digits = 3)

  sim_title <- simulation$get_name()
  tissue_title <- simulation$get_tissue_name()
  tissue_size <- paste(simulation$get_tissue_size(), collapse = " x ")

  #
  # cells %>%
  #   group_by(genotype, epistate) %>%
  #   summarise(n())

  color_map <- get_species_colors(simulation$get_species())

  ggplot2::ggplot(cells, ggplot2::aes(x = .data$position_x,
                                      y = .data$position_y,
                                      fill = .data$species)) +
    ggplot2::geom_hex(bins = num_of_bins) +
    my_theme() +
    ggplot2::labs(x = "Latitude", y = "Longitude",
                  fill = "Species",
                  caption = paste("Total number of cells",
                                  counts$counts %>% sum(),
                                  "(with", num_of_bins, "bins)"),
                  title = paste0(sim_title, " (t = ", time, ")"),
                  subtitle = paste("Tissue:", tissue_title, "[",
                                   tissue_size, "]")) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::xlim(0, simulation$get_tissue_size()[1]) +
    ggplot2::ylim(0, simulation$get_tissue_size()[2])
}



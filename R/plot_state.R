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

#' Plot the current state of the simulation.
#'
#' @description
#' A pie chart with population counts, split by species and epigentic state. It
#' also provides annotations for the simulation data.
#'
#' @param simulation A simulation.
#' @param color_map A named vector representing the simulation species color
#'   map (optional).
#' @return A ggplot plot.
#' @export
#'
#' @examples
#' sim <- SpatialSimulation()
#' sim$add_mutant(name = "A",
#'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
#'                growth_rates = c("+" = 0.2, "-" = 0.08),
#'                death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#'
#' plot_state(sim)
#'
#' # define a custom color map
#' color_map <- c("#B2DF8A", "#E31A1C")
#' names(color_map) <- c("A+", "A-")
#'
#' plot_state(sim, color_map=color_map)
plot_state <- function(simulation, color_map = NULL) {
  stopifnot(inherits(simulation, "Rcpp_SpatialSimulation"))

  counts <- simulation$get_counts() %>%
    dplyr::mutate(species = paste0(.data$mutant, .data$epistate))

  if (is.null(color_map)) {
    color_map <- get_species_colors(simulation$get_species())
  }

  ggplot2::ggplot(counts) +
    ggplot2::geom_bar(stat = "identity",
                      ggplot2::aes(x = "",
                                   y = .data$counts,
                                   fill = .data$species)) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(fill = "Species")
}
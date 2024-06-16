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

#' Plot the number of stochastic events in the simulation.
#'
#' @description
#' A pie chart with events split by type, mutant and epigentic state
#' where they occurred. It also provides annotations for the simulation
#' information.
#'
#' @param simulation A simulation.
#'
#' @return A ggplot plot.
#' @export
#' @examples
#' sim <- Simulation()
#' sim$add_mutant(name = "A",
#'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
#'                growth_rates = c("+" = 0.2, "-" = 0.08),
#'                death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#' plot_firings(sim)
plot_firings <- function(simulation) {
  stopifnot(inherits(simulation, "Rcpp_Simulation"))

  firings <- simulation$get_firings() %>%
    dplyr::mutate(species = paste0(.data$mutant, .data$epistate))

  time <- simulation$get_clock() %>% round(digits = 3)

  sim_title <- simulation$get_name()
  tissue_title <- simulation$get_tissue_name()
  tissue_size <- paste(simulation$get_tissue_size(), collapse = " x ")

  ggplot2::ggplot(firings) +
    ggplot2::geom_bar(stat = "identity",
                      ggplot2::aes(x = "", y = .data$fired,
                                   fill = .data$event)) +
    ggplot2::facet_grid(.data$mutant ~ .data$epistate) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(
      fill = "Event",
      x = "",
      y = "",
      title = paste0(sim_title, " (t = ", time, ")"),
      subtitle = paste(tissue_title, "of size", tissue_size),
      caption = paste("Total number of events", firings$fired %>% sum())
    ) +
    my_theme() +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::theme(legend.position = "bottom")
}
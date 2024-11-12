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
#' @param num_of_bins The number of bins (default: 100).
#' @param at_sample A sample name. This function will represent the tissue
#'   just before or just after the specified sampling (optional).
#' @param plot_sample_region A Boolean value. When `at_sample` is set,
#'   and `plot_sample_region` is set to be TRUE, this function also
#'   plots the sample region (default: TRUE).
#' @param color_map A named vector representing the simulation species color
#'   map (optional).
#' @param list_all_species A Boolean flag to show all species in
#'   the legend (default: FALSE).
#' @return An editable ggplot plot.
#' @export
#'
#' @examples
#' set.seed(0)
#' sim <- SpatialSimulation()
#' sim$add_mutant(name = "A",
#'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
#'                growth_rates = c("+" = 0.2, "-" = 0.08),
#'                death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_size("A-", 600)
#'
#' plot_tissue(sim)
#'
#' # define a custom color map
#' color_map <- c("#B2DF8A", "#E31A1C")
#' names(color_map) <- c("A+", "A-")
#'
#' plot_tissue(sim, color_map=color_map)
plot_tissue <- function(simulation, num_of_bins = 100,
                        at_sample = NULL,
                        plot_sample_region = TRUE,
                        color_map = NULL,
                        list_all_species = FALSE) {
  stopifnot(inherits(simulation, "Rcpp_SpatialSimulation"))

  # Get the cells
  if (is.null(at_sample)) {
    cells <- simulation$get_cells()
  } else {
    cells <- simulation$get_cells(at_sample)
  }

  cells <- cells %>% dplyr::as_tibble() %>%
    dplyr::mutate(species = paste0(.data$mutant, .data$epistate))

  if (is.null(color_map)) {
    color_map <- get_species_colors(simulation$get_species())
  }

  species_name_df <- simulation$get_species() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(species = paste0(.data$mutant, .data$epistate))

  cells$species <- factor(cells$species,
                          levels = unique(species_name_df$species))

  pl <- ggplot2::ggplot(cells, ggplot2::aes(x = .data$position_x,
                                            y = .data$position_y,
                                            fill = .data$species)) +
    ggplot2::geom_hex(bins = num_of_bins, show.legend = TRUE) +
    ggplot2::scale_fill_manual(values = color_map,
                               drop = !list_all_species) +
    my_theme() +
    ggplot2::labs(x = NULL, y = NULL,
                  fill = "Species") +
    # ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::xlim(-1, simulation$get_tissue_size()[1] + 1) +
    ggplot2::ylim(-1, simulation$get_tissue_size()[2] + 1)

  if (!is.null(at_sample) && plot_sample_region) {
    sample_info <- simulation$get_samples_info() %>%
      filter(.data$name == at_sample)

    pl <- pl + ggplot2::geom_rect(xmin = sample_info$xmin[[1]],
                                  xmax = sample_info$xmax[[1]],
                                  ymin = sample_info$ymin[[1]],
                                  ymax = sample_info$ymax[[1]],
                                  fill = NA, color = "black")
  }

  return(pl)
}



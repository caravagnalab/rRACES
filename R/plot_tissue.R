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
#' @param before_sample A sample name. When provided, this function represents
#'   the tissue as appeared just before the specified sampling. The parameters
#'   `before_sample` and `at_sample` are mutually exclusive (optional).
#' @param at_sample A sample name. When provided, this function represents
#'   the tissue as appeared when the specified sampling occurred. The
#'   parameters `before_sample` and `at_sample` are mutually exclusive
#'   (optional).
#' @param plot_next_sample_regions A Boolean value. When `before_sample` is
#'   set and `plot_next_sample_regions` is set to be TRUE, this function
#'   plots the regions of the samples collected at the same simulated time
#'   of the specified sample. When, instead, `at_sample` is set and
#'   `plot_next_sample_regions` is set to be TRUE, the function plots the
#'   regions of the samples collected at the same simulated time of the
#'   specified sample, but not before the specified sample (default: FALSE).
#' @param plot_sample_region A Boolean value. When either `at_sample` or
#'   `before_sample` are set and `plot_sample_region` is set to be TRUE,
#'   the function also plots the region of the specified sample
#'   (default: TRUE).
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
#' sim$run_up_to_size("A-", 60000)
#'
#' # collect 3 samples: "Sample_A", "Sample_B", and "Sample_C"
#' sim$sample_cells("Sample_A", c(425, 425), c(475, 475))
#' sim$sample_cells("Sample_B", c(525, 525), c(575, 575))
#' sim$sample_cells("Sample_C", c(425, 525), c(475, 575))
#'
#' # let the simulation evolve until the species "A-" account
#' # for 80000 cells
#' sim$run_up_to_size("A-", 80000)
#'
#' # plot the tissue in the current status
#' plot_tissue(sim)
#'
#' # plot the tissue as it was when "Sample_B" was collected
#' plot_tissue(sim, at_sample="Sample_B")
#'
#' # plot the tissue as it was when "Sample_B" was collected and
#' # highlight the regions of the samples collected at the same
#' # simulated time, but not before it, i.e., "Sample_B" and
#' # "Sample_C"
#' plot_tissue(sim, at_sample="Sample_B",
#'             plot_next_sample_regions = TRUE)
#'
#' # plot the tissue as it was just before sampling "Sample_B"
#' plot_tissue(sim, before_sample="Sample_B")
#'
#' # plot the tissue as it was just before sampling "Sample_B"
#' # and highlight the regions of the samples collected at the
#' # same simulated time, i.e., "Sample_A", "Sample_B", and
#' # "Sample_C"
#' plot_tissue(sim, before_sample="Sample_B",
#'             plot_next_sample_regions = TRUE)
#'
#' # define a custom color map
#' color_map <- c("#B2DF8A", "#E31A1C")
#' names(color_map) <- c("A+", "A-")
#'
#' plot_tissue(sim, color_map=color_map)
plot_tissue <- function(simulation, num_of_bins = 100,
                        before_sample = NULL,
                        at_sample = NULL,
                        plot_next_sample_regions = FALSE,
                        plot_sample_region = TRUE,
                        color_map = NULL,
                        list_all_species = FALSE) {
  stopifnot(inherits(simulation, "Rcpp_SpatialSimulation"))

  sample_info <- NULL

  # Get the cells
  if (is.null(at_sample)) {
    if (is.null(before_sample)) {
      cells <- simulation$get_cells()
    } else {
      samples <- simulation$get_samples_info()

      sample_info <- samples %>% filter(.data$name == before_sample)

      selected_samples <- samples %>%
        filter(.data$time == sample_info$time[[1]])

      first_sample_at_time <- samples %>%
        filter(.data$id == min(selected_samples$id))

      cells <- simulation$get_cells(first_sample_at_time$name[[1]])

      from_id <- first_sample_at_time$id[[1]]
    }
  } else {
    if (is.null(before_sample)) {
      samples <- simulation$get_samples_info()

      sample_info <- samples %>% filter(.data$name == at_sample)

      from_id <- sample_info$id[[1]]

      cells <- simulation$get_cells(at_sample)
    } else {
      stop(paste("The parameters `at_sample` and `before_sample`",
                 "are mutually exclusive."))
    }
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

  if (!is.null(sample_info)) {
    if (plot_sample_region) {
      pl <- pl + ggplot2::geom_rect(xmin = sample_info$xmin[[1]],
                                    xmax = sample_info$xmax[[1]],
                                    ymin = sample_info$ymin[[1]],
                                    ymax = sample_info$ymax[[1]],
                                    fill = NA, color = "black")
    }

    if (plot_next_sample_regions) {
      if (is.null(before_sample)) {
        selected_samples <- samples %>%
          filter(.data$time == sample_info$time[[1]],
                 .data$id >= from_id)
      }

      for (row_idx in seq_len(nrow(selected_samples))) {
        sample_info <- selected_samples[row_idx, ]
        pl <- pl + ggplot2::geom_rect(xmin = sample_info$xmin[[1]],
                                      xmax = sample_info$xmax[[1]],
                                      ymin = sample_info$ymin[[1]],
                                      ymax = sample_info$ymax[[1]],
                                      fill = NA, color = "black")
      }
    }
  }

  return(pl)
}



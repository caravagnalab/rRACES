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

#' Bounding box sampler
#'
#' @description
#' From a simulation object, this function samples bounding boxes of
#' a custom size until a number of cells of a certain species are not
#' sampled. If this never happens, after a fixed number of attempts
#' the best box is returned.
#'
#' @param simulation A simulation object.
#' @param which The species name.
#' @param n The desired number of cells from species `which`.
#' @param n_w Width of the box.
#' @param n_h Height of the box.
#' @param nattempts The maximum number of attempts.
#'
#' @return coordinates for a bounding box.
#' @export
#'
#' @examples
#' sim <- Simulation()
#' sim$add_mutant(name = "A", growth_rates = 0.08, death_rates = 0.01)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(60)
#' bbox <- bbox_sampler(sim, "A", n = 15, n_w = 5, n_h = 5)
#' sim$sample_cells("A", bbox$p, bbox$q)
#' plot_tissue(sim)
bbox_sampler <- function(simulation, which, n, n_w, n_h, nattempts = 100) {

  sp <- cli::make_spinner()

  # where are cells of species which
  locations <- simulation$get_cells() %>%
    dplyr::mutate(species = paste0(.data$mutant, .data$epistate)) %>%
    dplyr::filter(.data$species == which) %>%
    dplyr::select(.data$position_x, .data$position_y)

  bof_p <- bof_q <- NULL
  bof_counts <- -1

  repeat {
    sp$spin()

    location_id <- sample(seq_len(nrow(locations)), 1)

    p <- locations[location_id, ] %>% as.numeric()
    q <- p + c(n_w, n_h)

    nc <- simulation$get_cells(p, q) %>%
      dplyr::mutate(species = paste0(.data$mutant, .data$epistate)) %>%
      dplyr::filter(.data$species == which) %>%
      nrow()

    if (nc > bof_counts) {
      bof_p <- p
      bof_q <- q
      bof_counts <- nc
    }

    if (nc >= n || nattempts == 0) break
    nattempts <- nattempts - 1
  }

  sp$finish()

  return(list(p = bof_p, q = bof_q))
}

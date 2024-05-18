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

preprocess <- function(phylo_forest) {
  exposure <- phylo_forest$get_exposures()

  time_points <- exposure%>% dplyr::pull(time) %>% unique
  sbses <- exposure %>% dplyr::pull(signature) %>% unique
  for (t in time_points) {
    for (sbs in sbses) {
      if (nrow(exposure %>%
                 filter(.data$time == t,
                        .data$signature == sbs)) == 0) {
        exposure[nrow(exposure)+1,] <- c(as.numeric(t), sbs,
                                         as.numeric(0), "SBS")
      }
    }
  }

  last_time <- max(phylo_forest$get_samples_info()["time"])
  time_points <- sort(unique(append(time_points, last_time)))
  end_time <- apply(exposure, 1, function(x) {
    time_points[which(time_points == as.numeric(x[1])) + 1]
  })
  end_time[is.na(end_time)] <- as.numeric(last_time)

  exposure$end_time <- end_time

  exposure <- exposure[order(exposure$time, exposure$signature), ]

  exposure[, c(1, 3)] <- sapply(exposure[, c(1, 3)], as.numeric)

  return(exposure)
}

#' Plot a signatures exposures
#'
#' @description
#' Plots the signatures exposure changes along a phylogenetic
#' forest.
#'
#' @param phylogenetic_forest A phylogenetic forest.
#' @param linewidth The width of the lines in the plot.
#' @param emphatize_switches A Boolean flag to emphatize the
#'   exposure switches.
#'
#' @return An editable ggplot plot.
#' @examples
#' sim <- new(Simulation)
#'
#' sim$add_mutant(name = "A",
#'                growth_rates = 0.2,
#'                death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#'
#' sim$run_up_to_time(150)
#'
#' # sampling tissue
#'
#' n_w <- n_h <- 50
#' ncells <- 0.8 * n_w * n_h
#' bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
#'
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#'
#' forest <- sim$get_samples_forest()
#'
#' # placing mutations
#'
#' m_engine <- build_mutation_engine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name = "A",
#'                     passenger_rates = c(SNV = 8e-8))
#'
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8, ID3 = 1))
#' m_engine$add_exposure(time = 50,
#'                       c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5,
#'                         ID2 = 0.8, ID21 = 0.2))
#' phylo_forest <- m_engine$place_mutations(forest, 500)
#'
#' # plotting the phylogenetic forest
#' plot_exposure_timeline(phylo_forest)
#'
#' # plotting the phylogenetic forest emphatizing the exposure switches
#' plot_exposure_timeline(phylo_forest, emphatize_switches=TRUE)
#'
#' # deleting the mutation engine directory
#' unlink("demo", recursive=TRUE)
#' @export
plot_exposure_timeline <- function(phylogenetic_forest, linewidth = 0.8,
                                   emphatize_switches = FALSE) {
  stopifnot(inherits(phylogenetic_forest, "Rcpp_PhylogeneticForest"))

  exposure_df <- preprocess(phylogenetic_forest)

  signames <- exposure_df %>% dplyr::pull(signature) %>% unique
  all_colors <- rRACES:::get_signatures_colors()
  colors <- all_colors[signames]

  line_shift <- 0.005 * linewidth

  exposure_df <- exposure_df %>% dplyr::group_by(time,exposure) %>%
    dplyr::mutate(line_pos = .data$exposure + line_shift * (row_number() - 1))

  p <- ggplot2::ggplot(data = exposure_df,
                       ggplot2::aes(x = time, y = .data$line_pos,
                                    group = signature, color = signature)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$end_time,
                                       yend = .data$line_pos),
                          linewidth = linewidth) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = "Time Point", y = "Signature Exposure",
                  title = "Mutational Signature Exposure Timeline") +
    #scale_fill_manual(values=colors) +
    ggplot2::scale_colour_manual(name = "Signatures", values = colors) +
    ggplot2::theme_minimal()  # Apply a minimal theme

  if (emphatize_switches) {
    switching_times <- exposure_df["time"] %>% filter(time != 0)
    for (switching_time in switching_times) {
      p <- p +
        ggplot2::geom_vline(xintercept = switching_time, linetype = "dotted")
    }
  }

  return(p)
}

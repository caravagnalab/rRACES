#' Plot the current state of the simulation.
#'
#' @description
#' A pie chart with population counts, split by species and epigentic state. It
#' also provides annotations for the simulation data.
#'
#' @param simulation A simulation.
#'
#' @return A ggplot plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation, "plot_timeseries_test")
#' sim$history_delta <- 1 
#' sim$add_genotype(genotype = "A",
#'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
#'                  growth_rates = c("+" = 0.2, "-" = 0.08),
#'                  death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#'
#' # Set the frequency of storing
#' plot_timeseries(sim)
plot_timeseries <- function(simulation) {
  stopifnot(inherits(simulation, "Rcpp_Simulation"))

  counts <- simulation$get_count_history() %>%
    dplyr::mutate(species = paste0(.data$genotype, .data$epistate))
  time <- simulation$get_clock() %>% round(digits = 3)

  ncells <- simulation$get_counts() %>% dplyr::pull(counts) %>% sum()

  sim_title <- simulation$get_name()
  tissue_title <- simulation$get_tissue_name()
  tissue_size <- paste(simulation$get_tissue_size(), collapse = " x ")

  color_map <- get_species_colors(simulation$get_species())

  ggplot2::ggplot(counts) +
    ggplot2::geom_line(ggplot2::aes(x = time, y = .data$count,
                                    color = .data$species)) +
    ggplot2::geom_point(ggplot2::aes(x = time, y = .data$count,
                                     color = .data$species)) +
    ggplot2::labs(
      color = "Species",
      alpha = "Epistate", 
      title = paste0(sim_title, " (t = ", time, ")"),
      subtitle = paste("Tissue:", tissue_title, "[", tissue_size, "]"),
      caption = paste("Total number of cells", ncells)
    ) +
    my_theme() +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::theme(legend.position = "bottom")
}
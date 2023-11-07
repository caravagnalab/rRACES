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
#' sim <- new(Simulation, "plot_state_test")
#' sim$add_genotype(genotype = "A",
#'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
#'                  growth_rates = c("+" = 0.2, "-" = 0.08),
#'                  death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$add_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#' plot_state(sim)
plot_state <- function(simulation) {
  counts <- simulation$get_counts()
  time <- simulation$get_clock() %>% round(digits = 3)

  sim_title <- simulation$get_name()
  tissue_title <- simulation$get_tissue_name()
  tissue_size <- paste(simulation$get_tissue_size(), collapse = " x ")

  ggplot2::ggplot(counts) +
    ggplot2::geom_bar(stat = "identity",
                      ggplot2::aes(x = "", y = counts,
                                   fill = .data$genotype,
                                   alpha = .data$epistate)) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(
      fill = "Genotype",
      alpha = "Epistate", 
      title = paste0(sim_title, " (t = ", time, ")"),
      subtitle = paste("Tissue:", tissue_title, "[", tissue_size, "]"),
      caption = paste("Total number of cells", counts$counts %>% sum())
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::scale_fill_manual(values = get_species_colors(counts$genotype %>%
                                                             unique)) +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5))  +
    ggplot2::theme(legend.position = "bottom")
}
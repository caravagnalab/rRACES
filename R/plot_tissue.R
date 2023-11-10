#' Plot a tissue
#'
#' @description
#' Plots cells distribution over a tissue highlighting species by color.
#' To facilitate the plot and avoid excessive number of cells, for instance,
#' when a simulation deals with millions of cells, the plot draws a
#' hexagonal heatmap of 2D bins.
#'
#' @param simulation A simulation object.
#'
#' @return An editable ggplot plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation, "plot_tissue_test")
#' sim$add_genotype(name = "A",
#'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
#'                  growth_rates = c("+" = 0.2, "-" = 0.08),
#'                  death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#' plot_tissue(sim)
plot_tissue <- function(simulation) {
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

  ggplot2::ggplot(cells, ggplot2::aes(x = .data$position_x,
                                      y = .data$position_y,
                                      fill = .data$species)) +
    ggplot2::geom_hex(bins = 100) +
    my_theme() +
    ggplot2::labs(x = "Latitude", y = "Longitude",
                  fill = "Species",
                  caption = paste("Total number of cells",
                                  counts$counts %>% sum(), "(with 100 bins)"),
                  title = paste0(sim_title, " (t = ", time, ")"),
                  subtitle = paste("Tissue:", tissue_title, "[",
                                   tissue_size, "]")) +
    ggplot2::scale_fill_manual(values =
                                 get_species_colors(simulation)) +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::xlim(0, simulation$get_tissue_size()[1]) +
    ggplot2::ylim(0, simulation$get_tissue_size()[2])
}



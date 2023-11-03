#' Plot the number of stochastic events in the simulation.
#' 
#' @description
#' A piechart with events split by type, genotype and epigentic state where they
#' occurred. It also provides annotations for the simulation information.
#'
#' @param x A simulation.
#'
#' @return A ggplot plot.
#' @export
#' @examples
#' x = create_simulation()
#' x = add_species(x, "A")
#' x = set_initial_cell(x, "A", "+", c(50, 50))
#' x = run(x, list(time = 60))
#' plot_firings(x)
plot_firings = function(x)
{
  # t = x$get_simulated_clock()
  # counts = x$get_simulated_counts()
  # info = x$get_info()
  
  firings = x$simulation$get_firings() %>% dplyr::mutate(species = paste0(genotype, epistate))
  
  sim_title = x$name
  tissue_title = x$simulation$get_tissue_name()
  tissue_size = paste(x$simulation$get_tissue_size(), collapse = ' x ')
  
  ggplot2::ggplot(firings) +
    ggplot2::geom_bar(stat = 'identity', 
                      ggplot2::aes(x = "", y = fired, fill = event)) +
    ggplot2::facet_grid(genotype~epistate) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(
      fill = "Event", 
      x = "",
      y = "",
      title = paste0(sim_title, ' (t = 44)'),
      subtitle = paste(tissue_title, 'of size', tissue_size),
      caption = paste("Total number of events", firings$fired %>% sum())
    ) +
    # ggplot2::theme_void(base_size = 10) +
    my_theme() +
    ggplot2::scale_fill_brewer(palette = 'Dark2') +
    ggplot2::theme(legend.position = 'bottom')
}
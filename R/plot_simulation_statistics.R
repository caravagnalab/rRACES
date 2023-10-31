


#' Plot the number of stochastic events in the simulation.
#' 
#' @description
#' A piechart with events split by type, species and epigentic state where they
#' occurred. It also provides annotations for the simulation information.
#'
#' @param x A simulation.
#'
#' @return A ggplot plot.
#' @export
plot_firings = function(x)
{
  # t = x$get_simulated_clock()
  # counts = x$get_simulated_counts()
  # info = x$get_info()
  
  # Until the getters are fake
  firings = tibble::tibble(
    event = c("growth", "death", "switch+-", "switch-+", "growth", "death", "switch+-", "switch-+"),
    species = c("A", "A", "A", "A", "C", "C", "C", "C"),
    epistate = c("+", "+", "+", "-", "+", "+", "+", "-"),
    n = c(1000, 400, 3224, 334, 4223,2322, 2243,4224)
  )
  
  info = tibble::tibble(
    field = c("name", "tissue", "tissue_size"),
    value = c("My Simulation", "Liver", "100x100")
  )
  
  sim_title = info %>% filter(field == 'name') %>% pull(value)
  tissue_title = info %>% filter(field == 'tissue') %>% pull(value)
  tissue_size = info %>% filter(field == 'tissue_size') %>% pull(value)
  
  ggplot2::ggplot(firings) +
    ggplot2::geom_bar(stat = 'identity', 
                      ggplot2::aes(x = "", y = n, fill = event, alpha = epistate)) +
    ggplot2::facet_wrap(~species) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(
      fill = "Event", 
      alpha = 'Epigenetic', 
      title = paste0(sim_title, ' (t = 44)'),
      subtitle = paste(tissue_title, 'of size', tissue_size),
      caption = paste("Total number of events", firings$n %>% sum())
    ) +
    # ggplot2::theme_void(base_size = 10) +
    ggplot2::scale_fill_brewer(palette = 'Dark2') +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5))  +
    ggplot2::theme(legend.position = 'bottom')
}
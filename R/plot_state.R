#' Plot the current state of the simulation.
#' 
#' @description
#' A piechart with population counts, split by species and epigentic state. It 
#' also provides annotations for the simulation information.
#'
#' @param x A simulation.
#'
#' @return A ggplot plot.
#' @export
#'
#' @examples
plot_state = function(x)
{
  counts = x$simulation$get_counts()
  time = x$simulation$get_clock() %>% round(digits = 3)
  
  # # Until the getters are fake
  # counts = tibble::tibble(
  #   genotype = c("A", "A", "B", "B", "C", "C"),
  #   epistate = c("+", "-", "+", "-", "+", "-"),
  #   n = c(1000, 400, 3224, 334, 4223,2322)
  # )
  
  info = tibble::tibble(
    field = c("name", "tissue", "tissue_size"),
    value = c("My Simulation", "Liver", "100x100")
  )
  
  sim_title = info %>% filter(field == 'name') %>% pull(value)
  tissue_title = info %>% filter(field == 'tissue') %>% pull(value)
  tissue_size = info %>% filter(field == 'tissue_size') %>% pull(value)
  
  ggplot2::ggplot(counts) +
    ggplot2::geom_bar(stat = 'identity', 
                      ggplot2::aes(x = "", y = counts, fill = species, alpha = epistate)) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(
      fill = "Species", 
      alpha = 'Epigenetic', 
      title = paste0(sim_title, ' (t = ', time, ')'),
      subtitle = paste(tissue_title, 'of size', tissue_size),
      caption = paste("Total number of cells", counts$n %>% sum())
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::scale_fill_manual(values = get_species_colors(counts$species %>% unique)) +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5))  +
    ggplot2::theme(legend.position = 'bottom') 
}
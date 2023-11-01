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
plot_state = function(x)
{
  counts = x$simulation$get_counts()
  time = x$simulation$get_clock() %>% round(digits = 3)
  
  sim_title = x$name
  tissue_title = x$simulation$get_tissue_name()
  tissue_size = paste(x$simulation$get_tissue_size(), collapse = ' x ')
  
  ggplot2::ggplot(counts) +
    ggplot2::geom_bar(stat = 'identity', 
                      ggplot2::aes(x = "", y = counts, fill = genotype, alpha = epistate)) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(
      fill = "Genotype", 
      alpha = 'Epistate', 
      title = paste0(sim_title, ' (t = ', time, ')'),
      subtitle = paste("Tissue:", tissue_title, '[', tissue_size, ']'),
      caption = paste("Total number of cells", counts$n %>% sum())
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::scale_fill_manual(values = get_species_colors(counts$genotype %>% unique)) +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5))  +
    ggplot2::theme(legend.position = 'bottom') 
}
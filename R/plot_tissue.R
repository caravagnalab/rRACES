#' Plot a tissue
#' 
#' @description
#' Plots cells distribution over a tissue, higlightins species by colour. To facilitate
#' the plot and avoid eccessive number of cells (when there are for instance millions of
#' cells), the plot uses a 2D histogram.
#' 
#' @param x A simulation object.
#'
#' @return An editable ggplot plot.
#' @export
#'
#' @examples
#' x = create_simulation()
#' x = add_species(x, "A")
#' x = set_initial_cell(x, "A", "+", c(50, 50))
#' x = run(x, list(time = 60))
#' plot_tissue(x)
plot_tissue = function(x)
{
  stopifnot(inherits(x, 'rraces'))
  
  # Get the cells in the tissue at current simulation time
  cells = x$simulation$get_cells() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(species = paste0(genotype, epistate))
  
  counts = x$simulation$get_counts()
  
  time = x$simulation$get_clock() %>% round(digits = 3)
  
  sim_title = x$name
  tissue_title = x$simulation$get_tissue_name()
  tissue_size = paste(x$simulation$get_tissue_size(), collapse = ' x ')
  
  # 
  # cells %>%
  #   group_by(genotype, epistate) %>%
  #   summarise(n())
  
  ggplot2::ggplot(cells, ggplot2::aes(x = position_x, y = position_y, fill = genotype, alpha = epistate)) +
    ggplot2::geom_hex(bins = 100) +
    my_theme() +
    ggplot2::labs(x = "Latitude", y = 'Longitude', 
         fill = 'Genotype', alpha = 'Epistate',
         caption = paste("Total number of cells", counts$counts %>% sum(), "(with 100 bins)"),
         title = paste0(sim_title, ' (t = ', time, ')'),
         subtitle = paste("Tissue:", tissue_title, '[', tissue_size, ']')
    )+
    ggplot2::scale_fill_manual(values = get_species_colors(cells$genotype %>% unique)) +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5))  +
    ggplot2::theme(legend.position = 'bottom') +
    ggplot2::xlim(0, x$simulation$get_tissue_size()[1]) +
    ggplot2::ylim(0, x$simulation$get_tissue_size()[2])
}



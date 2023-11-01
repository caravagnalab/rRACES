#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
plot_tissue = function(x)
{
  # Get the cells in the tissue at current simulation time
  cells = x$simulation$get_cells() %>%
    as_tibble() %>%
    mutate(species = paste0(genotype, epistate))
  
  counts = x$simulation$get_counts()
  
  time = x$simulation$get_clock() %>% round(digits = 3)
  
  sim_title = x$name
  tissue_title = x$simulation$get_tissue_name()
  tissue_size = paste(x$simulation$get_tissue_size(), collapse = ' x ')
  
  # 
  # cells %>%
  #   group_by(genotype, epistate) %>%
  #   summarise(n())
  
  ggplot(cells, aes(x = position_x, y = position_y, fill = genotype, alpha = epistate)) +
    geom_hex(bins = 100) +
    my_theme() +
    labs(x = "Latitude", y = 'Longitude', 
         fill = 'Genotype', alpha = 'Epistate',
         caption = paste("Total number of cells", counts$counts %>% sum(), "(with 100 bins)"),
         title = paste0(sim_title, ' (t = ', time, ')'),
         subtitle = paste("Tissue:", tissue_title, '[', tissue_size, ']')
    )+
    ggplot2::scale_fill_manual(values = get_species_colors(cells$genotype %>% unique)) +
    ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5))  +
    ggplot2::theme(legend.position = 'bottom') +
    xlim(0, x$simulation$get_tissue_size()[1]) +
    ylim(0, x$simulation$get_tissue_size()[2])
}



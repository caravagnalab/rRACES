#' Annotate a plot of cell divisions.
#'
#' @description
#' It annotates a plot of cell divisions with information from sampling
#' times and MRCAs for all available samples
#'
#' @param tree_plot The output of `plot_forest`.
#' @param forest The original forest object from which the input to
#'        `plot_forest`
#' has been derived.
#' @param samples If `TRUE` it annotates samples.
#' @param MRCAs If `TRUE` it annotates MRCAs
#'
#' @return A `ggraph` tree plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation)
#' sim$add_mutant(name = "A", growth_rates = 0.08, death_rates = 0.01)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(60)
#' sim$sample_cells("MySample", c(500, 500), c(510, 510))
#' forest = sim$get_samples_forest()
#' forest$get_samples_info()
#' tree_plot = plot_forest(forest)
#' annotate_forest(tree_plot, forest)
annotate_forest <- function(tree_plot, forest, samples = TRUE, MRCAs = TRUE, 
                            exposures = FALSE, facet_signatures = TRUE) {
  
  samples_info <- forest$get_samples_info()
  
  # Sampling times
  if (samples) {
    
    max_Y <- max(tree_plot$data$y, na.rm = TRUE)
    
    tree_plot <- tree_plot +
      ggplot2::geom_hline(
        yintercept = max_Y - samples_info$time,
        color = "indianred3",
        linetype = "dashed",
        linewidth = .3
      )
  }
  
  # MRCAs
  if(MRCAs) {
    sample_names <- samples_info %>% dplyr::pull(.data$name)
    
    MRCAs_cells <- lapply(sample_names,
                          function(s){
                            forest$get_coalescent_cells(
                              forest$get_nodes() %>%
                                dplyr::filter(sample %in% s) %>%
                                dplyr::pull(.data$cell_id)
                            ) %>%
                              dplyr::mutate(sample = s)
                          }) %>%
      Reduce(f = dplyr::bind_rows) %>%
      dplyr::group_by(.data$cell_id) %>%
      dplyr::mutate(
        cell_id = paste(.data$cell_id)
      ) %>%
      dplyr::summarise(
        label = paste0("    ", .data$sample, collapse = "\n")
      )
    
    layout <- tree_plot$data %>%
      dplyr::select(.data$x, .data$y, .data$name) %>%
      dplyr::mutate(cell_id = paste(.data$name)) %>%
      dplyr::filter(.data$name %in% MRCAs_cells$cell_id) %>%
      dplyr::left_join(MRCAs_cells, by = "cell_id")
    
    
    tree_plot <-
      tree_plot +
      ggplot2::geom_point(
        data = layout,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = "purple3",
        size = 3,
        pch = 21
      ) +
      ggplot2::geom_text(
        data = layout,
        ggplot2::aes(x = .data$x, y = .data$y,
                     label = .data$label),
        color = "purple3",
        size = 3,
        hjust = 0,
        vjust = 1
      )
  }
  
  if(exposures){
    if(inherits(forest, "Rcpp_PhylogeneticForest")){
      max_Y <- max(tree_plot$data$y, na.rm = TRUE)
      # Get exposures table
      exposures <- forest$get_exposures()
      # Add exposures start and end times for each signature
      times <- exposures$time %>%  unique() %>%  sort()
      
      exposures <- exposures %>% 
        mutate(t_end = case_when(
          time == max(times) ~ Inf,
          TRUE ~ min(times[times > time]))
        ) %>% 
        mutate(signature = factor(signature, 
                                  levels = exposures %>% 
                                    dplyr::arrange(time) %>% 
                                    dplyr::pull(signature)))
      # Annotate exposures on tree
      tree_plot <- tree_plot +
        ggplot2::geom_rect(
          data = exposures,
          ggplot2::aes(
            xmin = -Inf,
            xmax = Inf,
            ymin = ifelse(is.infinite(t_end), 0, max_Y - t_end),
            ymax = max_Y - time,
            fill = signature,
            alpha = exposure
          )
        ) +
        ggplot2::scale_fill_manual(values = COSMIC_color_palette()) +
        ggplot2::scale_alpha_continuous(range = c(0.25, 0.75))+
        ggplot2::guides(fill = ggplot2::guide_legend(title = 'Signature'),
                        alpha = ggplot2::guide_legend(title = 'Exposure'))
      
      if(facet_signatures) tree_plot <- tree_plot + ggplot2::facet_wrap( ~ signature)
      # Push exposure rectangles to the back
      layers_new <- list(tree_plot$layers[[length(tree_plot$layers)]])
      layers_new <- c(layers_new, tree_plot$layers[1:length(tree_plot$layers)-1])
      
      tree_plot$layers <- layers_new
    }
    else cli::cli_alert_danger(text = "The input forest is not a PhylogeneticForest: cannot annotate signature exposures!")
  }
  
  tree_plot
}

#' @import Polychrome
COSMIC_color_palette <-  function(seed = 1999) {
  set.seed(seed)
  sigs = c('SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS10a', 'SBS10b', 'SBS11', 'SBS13', 'SBS14', 'SBS15', 'SBS17b', 'SBS18', 'SBS20', 'SBS21', 'SBS22a', 'SBS22b', 'SBS24', 'SBS25', 'SBS26', 'SBS28', 'SBS30', 'SBS31', 'SBS35', 'SBS36', 'SBS40a', 'SBS40b', 'SBS40c', 'SBS42', 'SBS87', 'SBS88', 'SBS90', 'SBS96', 'SBS97', 'SBS98', 'SBS99', 'SBS5', 'DBS1', 'DBS2', 'DBS3', 'DBS4', 'DBS5', 'DBS6', 'DBS7', 'DBS8', 'DBS9', 'DBS10', 'DBS11', 'DBS12', 'DBS13', 'DBS14', 'DBS15', 'DBS16', 'DBS17', 'DBS18', 'DBS19', 'DBS20', 'CN1', 'CN2', 'CN3', 'CN4', 'CN5', 'CN6', 'CN7', 'CN8', 'CN9', 'CN10', 'CN11', 'CN12', 'CN13', 'CN14', 'CN15', 'CN16', 'CN17', 'CN18', 'CN19', 'CN20', 'CN21', 'CN22', 'CN23', 'CN24', 'ID1', 'ID2', 'ID3', 'ID4', 'ID5', 'ID6', 'ID7', 'ID8', 'ID9', 'ID10', 'ID11', 'ID12', 'ID13', 'ID14', 'ID15', 'ID16', 'ID17', 'ID18', 'SBS8', 'SBS9', 'SBS10c', 'SBS10d', 'SBS12', 'SBS16', 'SBS17a', 'SBS19', 'SBS23', 'SBS27', 'SBS29', 'SBS32', 'SBS33', 'SBS34', 'SBS37', 'SBS38', 'SBS39', 'SBS41', 'SBS43', 'SBS44', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'SBS84', 'SBS85', 'SBS86', 'SBS89', 'SBS91', 'SBS92', 'SBS93', 'SBS94', 'SBS95') 
  return(Polychrome::createPalette(length(sigs), c("#6B8E23","#4169E1"), 
                                   M = 1000,
                                   target = "normal", range = c(15,80)) %>%
           setNames(sigs))
}

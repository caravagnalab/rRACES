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
#' sim$add_genotype(name = "A", growth_rates = 0.08, death_rates = 0.01)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(60)
#' sim$sample_cells("MySample", c(500, 500), c(510, 510))
#' forest = sim$get_samples_forest()
#' forest$get_samples_info()
#' tree_plot = plot_forest(forest)
#' annotate_forest(tree_plot, forest)
annotate_forest <- function(tree_plot, forest, samples = TRUE, MRCAs = TRUE) {
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
    sample_names <- samples_info %>% pull(.data$name)
    
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
  
  tree_plot
}
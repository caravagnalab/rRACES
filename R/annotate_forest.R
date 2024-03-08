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
#' @param exposures If `TRUE` it annotates exposures to mutational signatures
#' @param facet_signatures If `TRUE` and if `exposures` is `TRUE` it creates a faceted forest plot where the exposure to each signature is annotated on a separated plot 
#' @param drivers If `TRUE` it annotates drivers on the node they originated.
#' @param add_driver_label If `TRUE` and if `drivers` is `TRUE` it annotates the driver name.
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
                            exposures = FALSE, facet_signatures = TRUE, 
                            drivers = TRUE, add_driver_label = TRUE) {
  
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
  
  if(exposures) {
    if(inherits(forest, "Rcpp_PhylogeneticForest")) {
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
        ggplot2::scale_fill_manual(values = get_signatures_colors()) +
        ggplot2::scale_alpha_continuous(range = c(0.25, 0.75),
                                        breaks = sort(unique(exposures$exposure))) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Signature"),
                        alpha = ggplot2::guide_legend(title = "Exposure"))
      
      if(facet_signatures) tree_plot <- tree_plot + ggplot2::facet_wrap( ~ signature)
      # Push exposure rectangles to the back
      layers_new <- list(tree_plot$layers[[length(tree_plot$layers)]])
      layers_new <- c(layers_new, tree_plot$layers[1:length(tree_plot$layers)-1])
      
      tree_plot$layers <- layers_new
    }
    else cli::cli_alert_danger(text = "The input forest is not a PhylogeneticForest: cannot annotate signature exposures!")
  }
  
  if(drivers) {
    if(inherits(forest, "Rcpp_PhylogeneticForest")) {
      drivers_SNVs = drivers_CNAs = data.frame()
      try(expr = {
        drivers_SNVs <- forest$get_sampled_cell_SNVs() %>% 
          dplyr::filter(class == "driver") %>% 
          dplyr::mutate(driver_id = paste0(chromosome, ":", chr_pos, ":", ref, ">", alt),
                        driver_type = "SNV") %>% 
          dplyr::select(cell_id, driver_id, driver_type)
      })
      try(expr = {
        drivers_CNAs <- forest$get_sampled_cell_CNAs() %>% 
          dplyr::mutate(driver_id = paste0(chromosome, ":", begin, "-", end, ":", allele),
                        driver_type = "CNA") %>% 
          dplyr::select(cell_id, driver_id, driver_type)
      })
      
      drivers <- dplyr::bind_rows(drivers_SNVs, drivers_CNAs)
      
      drivers_start_nodes <- lapply(unique(drivers$driver_id), function(d) {
        nodes_with_driver = drivers %>% dplyr::filter(driver_id==d) %>% dplyr::pull(cell_id)
        d_type = drivers %>% dplyr::filter(driver_id==d) %>% dplyr::pull(driver_type) %>% unique()
        
        forest$get_coalescent_cells(nodes_with_driver) %>% 
          dplyr::mutate(driver_id=d, driver_type=d_type)
      }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(cell_id = as.character(cell_id)) %>% 
        dplyr::group_by(cell_id) %>% 
        dplyr::summarise(driver_id = paste0(driver_id, collapse = "\n"))
      
      layout <- tree_plot$data %>%
        dplyr::select(x, y, name) %>%
        dplyr::mutate(cell_id = paste(name), has_driver=TRUE) %>%
        dplyr::filter(name %in% drivers_start_nodes$cell_id) %>%
        dplyr::left_join(drivers_start_nodes, by = "cell_id")
      
      tree_plot <- tree_plot +
        ggplot2::geom_point(
          data = layout,
          ggplot2::aes(x = .data$x, y = .data$y),
          fill = "#d68910",
          color = "#FF000000",
          size = 2,
          pch = 21
        )
      
      if(add_driver_label) {
        nudge_x = (max(tree_plot$data$x) - min(tree_plot$data$x)) * .15
        tree_plot <- tree_plot +
          ggrepel::geom_label_repel(
            data = layout,
            ggplot2::aes(x = .data$x, y = .data$y,
                         label = .data$driver_id),
            color = "#d68910",
            size = 2.5,
            hjust = 0,
            nudge_x = -nudge_x,
            direction = "x"
          )
      } 
    } 
    else cli::cli_alert_danger(text = "The input forest is not a PhylogeneticForest: cannot annotate signature exposures!")
  }
  
  tree_plot
}

#' Annotate a plot of cell divisions.
#'
#' @description
#' It annotates a plot of cell divisions where branches containing relevant biological events are
#' colored
#'
#' @param forest The original forest object
#' has been derived.
#' @param labels Table annotating the sticks (output of get_events_table).
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
#' labels = get_events_table(forest)
#' plot_sticks(forest, labels)

plot_sticks = function(forest, labels){
  stopifnot(inherits(forest, "Rcpp_SamplesForest"))
  nodes <- forest$get_nodes()
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  }
  else {
    
    forest_data <- forest$get_nodes() %>% dplyr::as_tibble()  %>%
      dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
      dplyr::select(.data$from,
                    .data$to,
                    .data$mutant,
                    .data$epistate,
                    .data$sample,
                    .data$birth_time) %>%
      dplyr::mutate(
        from = ifelse(is.na(.data$from), "WT",
                      .data$from),
        species = paste0(.data$mutant,
                         .data$epistate),
        sample = ifelse(is.na(.data$sample),
                        "N/A", .data$sample),
        highlight = FALSE
      ) %>%
      full_join(labels %>% 
                  dplyr::rename(to = cell_id) %>% 
                  dplyr::select(c("to","obs")), by = "to") %>% 
       mutate(obs = ifelse(to == 0,"pre-neoplastic",obs))
    
    edges <- forest_data %>% dplyr::select("from", "to", 
                                           "highlight")
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)
    graph <- graph %>% tidygraph::activate("nodes") %>% 
      dplyr::left_join(forest_data %>% dplyr::rename(name = .data$to) %>% 
                         dplyr::mutate(name = as.character(.data$name)), 
                       by = "name")
    layout <- ggraph::create_layout(graph, layout = "tree", 
                                    root = "WT") 
    max_Y <- max(layout$birth_time, na.rm = TRUE)
    layout$reversed_btime <- max_Y - layout$birth_time
    layout$y <- layout$reversed_btime
    # species_colors <- get_species_colors(forest$get_species_info())
    cls = ggsci::pal_simpsons()(forest_data %>% 
                                  filter(obs != "pre-neoplastic") %>% pull(obs) %>% unique() %>% length())
    names(cls) = forest_data %>% filter(obs != "pre-neoplastic")  %>% arrange(desc(to)) %>% 
      pull(obs) %>% unique()
    cls["Subclonal"] = "gainsboro"
    cls["pre-neoplastic"] = "salmon"
    ncells <- graph %>% tidygraph::activate(nodes) %>% dplyr::pull(.data$name) %>% 
      length()
    ncells_sampled <- graph %>% tidygraph::activate(nodes) %>% 
      dplyr::filter(sample != "N/A") %>% dplyr::pull(.data$name) %>% 
      length()
    nsamples <- forest$get_samples_info() %>% nrow()
    labels_every <- max_Y/10
    point_size <- c(0,rep(1,length(forest_data$obs[forest_data$obs != "Subclonal"] %>% unique())))
    names(point_size) <- c("Subclonal", forest_data$obs[forest_data$obs != "Subclonal"] %>% unique())
    
    graph_plot = ggraph::ggraph(layout, "tree") + 
      ggraph::geom_edge_link(edge_width = 0.5,
                             ggplot2::aes(edge_color = ifelse(highlight, "indianred3","gainsboro")))
    # graph_plot + ggraph::geom_node_point(ggplot2::aes(color = .data$species,
    graph_plot + ggraph::geom_node_point(ggplot2::aes(
      color = .data$obs,
      #, shape = ifelse(is.na(.data$sample), "N/A", .data$sample),
      size = .data$obs
    )) +
      ggplot2::scale_color_manual(values = cls, na.translate=FALSE)  + 
      # ggplot2::scale_color_manual(values = species_colors) +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        title = "Phylogenetic length",
        subtitle = paste0(
          ncells,
          " cells (",
          ncells_sampled,
          " sampled from ",
          nsamples,
          " samples)"
        ),
        shape = "Sample",
        x = NULL,
        y = "Birth time"
      ) + ggplot2::guides(
        size = "none",
        shape = ggplot2::guide_legend("Sample"),
        fill = ggplot2::guide_legend("Species")
      ) +
      ggplot2::scale_size_manual(values = point_size) +
      ggplot2::scale_y_continuous(   labels = seq(0, max_Y,
                                                  labels_every) %>% round %>% rev,
                                     breaks = seq(0,
                                                  max_Y, labels_every) %>% round
      ) + ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
}

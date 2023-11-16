#' Plot a forest of cell divisions.
#' 
#' @description
#' Plot a forest of trees for cells passed in input, after a call to method
#' `get_nodes` from class `SamplesForest`. This plot is carried out using
#' `ggraph` and for simplicity of visualisation the forest is plot as a
#' set of trees connected to a generic wildtype cell.
#'
#' @param forest_nodes Result from`get_nodes` from class `SamplesForest`.
#'
#' @return A `ggraph` tree plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation, "plot_sampled_cells_test")
#' sim$add_genotype(name = "A", growth_rates = 0.08, death_rates = 0.01)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(60)
#' sim$sample_cells("MySample", c(500, 500), c(510, 510))
#' forest = sim$get_samples_forest()
#' forest$get_samples_info()
#' plot_sampled_cells(forest$get_nodes())
plot_sampled_cells = function(forest_nodes)
{
  forest_data = forest_nodes %>% 
    dplyr::as_tibble() %>% 
    dplyr::rename(
      from = ancestor,
      to = cell_id
    ) %>% 
    dplyr::select(from, to, genotype, epistate, sample) %>% 
    dplyr::mutate(
      from = ifelse(is.na(from), "WT", from),
      species = paste0(genotype, epistate),
      sample = ifelse(is.na(sample), "N/A", sample)
    )
  
  edges = forest_data %>% 
    dplyr::select(from, to)
  
  # Create a tidygraph object
  graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)
  
  graph = graph %>% 
    tidygraph::activate(nodes) %>% 
    dplyr::left_join(
      forest_data %>% 
        dplyr::rename(name = to) %>% 
        dplyr::mutate(name = as.character(name)),
      by = 'name'
    )
  
  # Define the layout with the tree rooted and the root on top
  layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")
  
  species_colors = rRACES:::get_species_colors(sim)
  
  ncells = graph %>%
    tidygraph::activate(nodes) %>% 
    dplyr::pull(name) %>% 
    length
  
  ncells_sampled = graph %>%
    tidygraph::activate(nodes) %>% 
    dplyr::filter(sample != "N/A") %>% 
    dplyr::pull(name) %>% 
    length
  
  nsamples = forest$get_samples_info() %>% nrow()
  
  labels_every = (layout$y %>% max)/10
  
  point_size = c(.5, rep(1, nsamples))
  names(point_size) = c("N/A", forest$get_samples_info() %>% pull(name))
  
  # sample_times = forest$get_samples_info() %>% 
  #   pull(time) %>% 
  #   unique
  # 
  # point_shape = c(.5, rep(1, nsamples))
  # names(point_size) = c("N/A", forest$get_samples_info() %>% pull(name))
  
  # Plot the graph with color-coded nodes
  # figure = 
  ggraph::ggraph(layout, "tree") +
    ggraph::geom_edge_link(edge_width = .1) +
    ggraph::geom_node_point(
      ggplot2::aes(
        color = species, 
        shape = ifelse(is.na(sample), "N/A", sample),
        size = sample
      )) +
    ggplot2::scale_color_manual(values = species_colors) +
    # scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    # theme_void() +
    ggplot2::theme_minimal() + 
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = "Cell divisons",
      subtitle = paste0(ncells, " cells (", ncells_sampled, ' sampled from ', nsamples, " samples)"),
      shape = "Sample",
      x = NULL,
      y = "Cell division"
    ) +
    ggplot2::guides(size = 'none', fill = "Species") +
    ggplot2::scale_size_manual(
      values = point_size
    ) +
    ggplot2::scale_y_continuous(labels = 
                         seq(0, layout$y %>% max, labels_every) %>% round %>% rev,
                       breaks = 
                         seq(0, layout$y %>% max, labels_every) %>% round
    ) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
}

# ggsave(
#   "~/Downloads/tree.pdf",
#   figure,
#   width = 30,
#   height = 20
# )
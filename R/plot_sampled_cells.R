#' Plotcell divisions.
#' 
#' @description
#' Plot a forest of trees for cells passed in input, as those that can be 
#' obtained with `get_nodes` from class `SamplesForest`. This plot is carried out using
#' `ggraph`.
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
    dplyr::select(from, to, genotype, epistate, sample, birth_time) %>% 
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
  
  max_Y = max(layout$birth_time, na.rm = TRUE)
  layout$reversed_btime = max_Y - layout$birth_time
  
  layout$y = layout$reversed_btime

  species_colors = get_species_colors(sim)
  
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
  
  labels_every = max_Y/10
  
  point_size = c(.5, rep(1, nsamples))
  names(point_size) = c("N/A", forest$get_samples_info() %>% pull(name))
  
  # Plot the graph with color-coded nodes
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
      title = "Birth time",
      subtitle = paste0(ncells, " cells (", ncells_sampled, ' sampled from ', nsamples, " samples)"),
      shape = "Sample",
      x = NULL,
      y = "Cell division"
    ) +
    ggplot2::guides(size = 'none', 
                    shape = ggplot2::guide_legend("Sample"),
                    fill = ggplot2::guide_legend("Species")) +
    ggplot2::scale_size_manual(
      values = point_size
    ) +
    ggplot2::scale_y_continuous(labels = 
                         seq(0, max_Y, labels_every) %>% round %>% rev,
                       breaks = 
                         seq(0, max_Y, labels_every) %>% round
    ) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
}

#' Annotate a plot of cell divisions.
#' 
#' @description
#' It annotates a plot of cell divisions with information from sampling
#' times and MRCAs for all available samples
#'
#' @param tree_plot The output of `plot_sampled_cells`.
#' @param forest The original forest object from which the input to `plot_sampled_cells`
#' has been derived.
#' @param samples If `TRUE` it annotates samples.
#' @param MRCAs If `TRUE` it annotates MRCAs
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
#' tree_plot = plot_sampled_cells(forest$get_nodes())
#' annotate_forest(tree_plot, forest)
annotate_forest = function(tree_plot, forest, samples = TRUE, MRCAs = TRUE)
{
  forest_nodes = forest$get_nodes()
  # Sampling times
  if(samples)
  {
    samples_table = forest$get_samples_info()
    
    max_Y = max(tree_plot$data$y, na.rm = TRUE)
    
    tree_plot = tree_plot +
      ggplot2::geom_hline(
        yintercept = max_Y - samples_table$time,
        color = 'indianred3',
        linetype = 'dashed',
        linewidth = .3
      )
  }
  
  # MRCAs
  if(MRCAs)
  {
    samples_list = forest$get_samples_info() %>% pull(name)
    
    MRCAs_cells = lapply(
      samples_list, 
      function(s){
        forest$get_coalescent_cells(
          forest_nodes %>% 
            dplyr::filter(sample %in% s) %>% 
            dplyr::pull(cell_id)
        ) %>% 
          dplyr::mutate(sample = s)
      }) %>% 
      Reduce(f = dplyr::bind_rows) %>% 
      dplyr::group_by(cell_id) %>% 
      dplyr::mutate(
        cell_id = paste(cell_id)
      ) %>% 
      dplyr::summarise(
        label = paste0("    ", sample, collapse = '\n')
      )
    
    layout = tree_plot$data %>% 
      dplyr::select(x, y, name) %>% 
      dplyr::mutate(cell_id = paste(name)) %>%  
      dplyr::filter(name %in% MRCAs_cells$cell_id) %>% 
      dplyr::left_join(MRCAs_cells, by = 'cell_id')
    
    tree_plot = 
      tree_plot +
      ggplot2::geom_point(
        data = layout,
        ggplot2::aes(x = x, y = y),
        color = 'purple3',
        size = 3,
        pch = 21
      ) +
      ggplot2::geom_text(
        data = layout,
        ggplot2::aes(x = x, y = y, label = label),
        color = 'purple3',
        size = 3,
        hjust = 0,
        vjust = 1
      ) 
  }
  
  tree_plot
}

# ggsave(
#   "~/Downloads/tree.pdf",
#   figure,
#   width = 30,
#   height = 20
# )
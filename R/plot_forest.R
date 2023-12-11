## This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
## Copyright (C) 2023 - Giulio Caravagna <gcaravagna@units.it>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' Plot samples cell division forest.
#'
#' @description
#' Plot the samples cell division forest. This plot is carried out using
#' `ggraph` and for simplicity of visualisation the forest is plot as a
#' set of trees connected to a generic wildtype cell.
#'
#' @param forest The samples forest to be plot.
#' @param highlight_sample If a sample name, the path from root to the sampled
#' cells in the sample is highlighted. If `NULL` (default), nothing is highlighted.
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
#' forest <- sim$get_samples_forest()
#' plot_forest(forest)
plot_forest <- function(forest, highlight_sample = NULL) {
  stopifnot(inherits(forest, "Rcpp_SamplesForest"))

  nodes <- forest$get_nodes()

  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    forest_data <- forest$get_nodes() %>%
      dplyr::as_tibble() %>%
      dplyr::rename(
        from = .data$ancestor,
        to = .data$cell_id
      ) %>%
      dplyr::select(.data$from, .data$to, .data$mutant,
                    .data$epistate, .data$sample, .data$birth_time) %>%
      dplyr::mutate(
        from = ifelse(is.na(.data$from), "WT", .data$from),
        species = paste0(.data$mutant, .data$epistate),
        sample = ifelse(is.na(.data$sample), "N/A", .data$sample),
        highlight = FALSE
      )
    
    # Highlight is optional
    if(!is.null(highlight_sample))
    {
      highlight = paths_to_sample(forest_data, highlight_sample)
      forest_data$highlight = forest_data$to %in% highlight$to
    }

    edges <- forest_data %>%
      dplyr::select("from", "to", "highlight")

    # Create a tidygraph object
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)

    graph <- graph %>%
      tidygraph::activate("nodes") %>%
      dplyr::left_join(
        forest_data %>%
          dplyr::rename(name = .data$to) %>%
          dplyr::mutate(name = as.character(.data$name)),
        by = "name"
      )

    # Define the layout with the tree rooted and the root on top
    layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")

    max_Y <- max(layout$birth_time, na.rm = TRUE)
    layout$reversed_btime <- max_Y - layout$birth_time

    layout$y <- layout$reversed_btime

    species_colors <- get_species_colors(forest$get_species_info())

    ncells <- graph %>%
      tidygraph::activate(nodes) %>%
      dplyr::pull(.data$name) %>%
      length()

    ncells_sampled <- graph %>%
      tidygraph::activate(nodes) %>%
      dplyr::filter(sample != "N/A") %>%
      dplyr::pull(.data$name) %>%
      length()

    nsamples <- forest$get_samples_info() %>% nrow()

    labels_every <- max_Y/10

    point_size <- c(.5, rep(1, nsamples))
    names(point_size) <- c("N/A", forest$get_samples_info() %>%
                            pull(.data$name))
    
    # highlight_edges = paths_to_sample(layout, "S_2_2")
    
    # layout$highlight = layout$name %in% highlight_edges$name

    # Plot the graph
    graph_plot = ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = .1, 
                             ggplot2::aes(edge_color = ifelse(highlight, 
                                                              "indianred3", 
                                                              "black")))

    graph_plot +
      ggraph::geom_node_point(ggplot2::aes(color = .data$species,
                                           shape = ifelse(is.na(.data$sample),
                                                          "N/A",
                                                          .data$sample),
                                           size = .data$sample)) +
      ggplot2::scale_color_manual(values = species_colors) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(title = "Birth time",
        subtitle = paste0(ncells, " cells (", ncells_sampled,
                          " sampled from ", nsamples, " samples)"),
        shape = "Sample",
        x = NULL,
        y = "Cell division"
      ) +
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend("Sample"),
                      fill = ggplot2::guide_legend("Species")) +
      ggplot2::scale_size_manual(
        values = point_size
      ) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y,
                                               labels_every) %>%
          round %>% rev,
        breaks = seq(0, max_Y, labels_every) %>% round
      ) +
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
}


paths_to_sample = function(forest_data, sample){
  
  to = forest_data %>% dplyr::filter(sample == !!sample) %>% pull(to)
  
  Q = NULL
  
  while(length(to) > 0)
  {
    # New element to inspect
    to_head = to[1]
    to = to[-1]

    # Forward star 
    to_add = forest_data %>%  dplyr::filter(to == to_head)
    Q = dplyr::bind_rows(Q, to_add) %>% dplyr::distinct()
    
    # Recursion
    to = c(to, to_add %>% dplyr::pull(.data$from)) %>% unique()
  }
  
  Q
}


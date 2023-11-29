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
#' forest <- sim$get_samples_forest()
#' plot_forest(forest)
plot_forest <- function(forest) {
  stopifnot(inherits(forest, "Rcpp_SamplesForest"))

  nodes <- forest$get_nodes()

  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot())
  } else {
    forest_data <- forest$get_nodes() %>%
      dplyr::as_tibble() %>%
      dplyr::rename(
        from = .data$ancestor,
        to = .data$cell_id
      ) %>%
      dplyr::select(.data$from, .data$to, .data$genotype,
                    .data$epistate, .data$sample, .data$birth_time) %>%
      dplyr::mutate(
        from = ifelse(is.na(.data$from), "WT", .data$from),
        species = paste0(.data$genotype, .data$epistate),
        sample = ifelse(is.na(.data$sample), "N/A", .data$sample)
      )

    edges <- forest_data %>%
      dplyr::select(.data$from, .data$to)

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

    # Plot the graph with color-coded nodes
    ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = .1) +
      ggraph::geom_node_point(ggplot2::aes(color = .data$species,
                                           shape = ifelse(is.na(.data$sample),
                                                          "N/A",
                                                          .data$sample),
                                           size = .data$sample)) +
      ggplot2::scale_color_manual(values = species_colors) +
      # scale_shape_manual(values = c("+" = 16, "-" = 17)) +
      # theme_void() +
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

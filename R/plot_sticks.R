#' Annotate a plot of cell divisions.
#'
#' @description
#' It annotates a plot of cell divisions where branches containing relevant
#'   biological events are colored
#'
#' @param forest The original forest object has been derived.
#' @param labels A data frame annotating the sticks (it can be the output of
#'   `get_relevant_branches()`).
#' @param cls A custom list of colors for any stick. If NULL a deafult palette
#'   is chosen.
#'
#' @return A `ggraph` tree plot.
#' @export
#'
#' @examples
#' sim <- SpatialSimulation()
#'
#' sim$add_mutant(name = "A",
#'                growth_rates = 1,
#'                death_rates = 0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_size("A",1e4)
#'
#' sim$add_mutant(name = "B",
#'                growth_rates = 3.5,
#'                death_rates = 0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#'
#' sim$run_up_to_size("B",1e4)
#'
#' bbox <- sim$search_sample(c("A" = 100,"B" = 100), 50, 50)
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#' forest = sim$get_samples_forest()
#'
#' labels = get_relevant_branches(forest)
#' plot_sticks(forest, labels)

plot_sticks = function(forest, labels, cls = NULL) {
  stopifnot(inherits(forest, "Rcpp_SamplesForest"))
  nodes <- forest$get_nodes()
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    labels <- labels %>%
      dplyr::mutate(
        cell_id = as.character(.data$cell_id),
        ancestor = ifelse(is.na(.data$ancestor), "WT",
                          as.character(.data$ancestor))
      ) %>%
      dplyr::add_row(cell_id="WT", ancestor="WT", mutant=NA, epistate="",
                     sample=NA, birth_time=0, label="Truncal")

    forest_data <- forest$get_nodes()

    forest_data[nrow(forest_data) + 1, ] <- c(NA, NA, NA, NA, NA, 0)

    forest_data <- forest_data %>% dplyr::as_tibble()  %>%
      dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
      dplyr::select(.data$from,
                    .data$to,
                    .data$mutant,
                    .data$epistate,
                    .data$sample,
                    .data$birth_time) %>%
      dplyr::mutate(
        from = ifelse(is.na(.data$from), "WT", .data$from),
        to = ifelse(is.na(.data$to), "WT", as.character(.data$to)),
        species = paste0(.data$mutant, .data$epistate),
        sample = ifelse(is.na(.data$sample), "N/A", .data$sample),
        highlight = FALSE
      ) %>%
      dplyr::full_join(labels %>%
                  dplyr::rename(to = cell_id) %>%
                  dplyr::select(c("to", "label")),
                by = "to") %>%
      dplyr::mutate(label = ifelse(is.na(label), "Subclonal", label))

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

    if (is.null(cls)) {
      cls <- get_colors_for(forest_data %>% dplyr::arrange(desc(to)) %>%
                                dplyr::pull(label) %>% unique())
      cls["Subclonal"] <- "gainsboro"
    }
    ncells <- graph %>% tidygraph::activate(nodes) %>% dplyr::pull(.data$name) %>%
      length()
    ncells_sampled <- graph %>% tidygraph::activate(nodes) %>%
      dplyr::filter(sample != "N/A") %>% dplyr::pull(.data$name) %>%
      length()
    nsamples <- forest$get_samples_info() %>% nrow()
    labels_every <- max_Y/10
    not_subclonal <- forest_data$label[forest_data$label != "Subclonal"] %>%
        unique()
    point_size <- c(0,rep(1,length(not_subclonal)))
    names(point_size) <- c("Subclonal", not_subclonal)

    ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = 0.5,
                             ggplot2::aes(edge_color = ifelse(highlight,
                                                              "indianred3",
                                                              "gainsboro"))) +

      ggraph::geom_node_point(ggplot2::aes(
        color = .data$label,
        size = .data$label
      )) +
      ggplot2::scale_color_manual(values = cls, na.translate=FALSE)  +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        title = "Phylogenetic tree",
        subtitle = paste0(ncells, " cells (", ncells_sampled,
                          " sampled from ", nsamples, " samples)"
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
      ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>%
                                    round %>% rev,
                                  breaks = seq(0, max_Y, labels_every) %>%
                                    round
      ) + ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
}

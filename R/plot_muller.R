collapse_loops <- function(df_edges) {
  queue <- c("Wild-type")
  Q <- collapsed_graph <- NULL
  while (length(queue) > 0) {
    parent <- queue[1]
    Q <- c(Q, parent)

    w_edg <- df_edges %>% dplyr::filter(.data$Parent == parent,
                                        !(.data$Identity %in% Q))

    collapsed_graph <- dplyr::bind_rows(collapsed_graph, w_edg)

    queue <- queue[-1]
    queue <- c(queue, w_edg$Identity)
  }

  collapsed_graph
}

#' Plot a muller plot
#'
#' @description
#' Plots a muller plot of the simulation, using the store time-series data 
#' and package ggmuller.
#'
#' @param simulation A simulation object.
#'
#' @return An editable ggplot plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation, "plot_muller_test")
#' sim$add_genotype(name = "A",
#'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
#'                  growth_rates = c("+" = 0.2, "-" = 0.08),
#'                  death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$history_delta = 1
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#' plot_muller(sim)
plot_muller <- function(simulation) {
  stopifnot(inherits(simulation, "Rcpp_Simulation"))

  # Tumour DF
  df_populations <- simulation$get_count_history() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Identity = paste0(.data$genotype, .data$epistate)) %>%
    dplyr::rename(Generation = .data$time, Population = .data$count) %>%
    dplyr::select(.data$Generation, .data$Identity, .data$Population)

  # Tumour edges
  df_edges <- simulation$get_lineage_graph() %>%
    # dplyr::mutate(lblr = paste(ancestor, progeny),
    # lbrl = paste(progeny, ancestor)) %>%
    dplyr::distinct(.data$ancestor, .data$progeny) %>%
    dplyr::rename(Parent = .data$ancestor, Identity = .data$progeny) %>%
    dplyr::select(.data$Parent, .data$Identity)

  # Remove loops
  df_edges <- collapse_loops(df_edges)

  max_tumour_size <- df_populations %>%
    dplyr::group_by(.data$Generation) %>%
    dplyr::summarise(Population = sum(.data$Population)) %>%
    dplyr::pull(.data$Population) %>%
    max()

  max_tumour_size <- max_tumour_size * 1.05

  # Wild-type dynamics
  wt_dynamics <- df_populations %>%
    dplyr::group_by(.data$Generation) %>%
    dplyr::summarise(Population = sum(.data$Population)) %>%
    dplyr::mutate(Identity = "Wild-type",
                  Population = max_tumour_size - .data$Population)

  t_wt_dynamics <- dplyr::bind_rows(
    wt_dynamics,
    df_populations
  )

  muller_df <- ggmuller::get_Muller_df(df_edges, t_wt_dynamics)

  color_map <- get_species_colors(simulation$get_species())

  ggmuller::Muller_pop_plot(muller_df, add_legend = TRUE) +
    my_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend("Species")) +
    ggplot2::scale_fill_manual(values = c(`Wild-type` = "gainsboro",
                                          color_map)
    )
}

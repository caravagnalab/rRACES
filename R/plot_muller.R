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
#' sim <- new(Simulation, "plot_tissue_test")
#' sim$add_genotype(name = "A",
#'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
#'                  growth_rates = c("+" = 0.2, "-" = 0.08),
#'                  death_rates = c("+" = 0.1, "-" = 0.01))
#' sim$history_delta = 1
#' sim$place_cell("A+", 500, 500)
#' sim$run_up_to_time(60)
#' plot_muller(sim)
plot_muller = function(simulation)
{
  stopifnot(inherits(simulation, "Rcpp_Simulation"))
  
  # Tumour DF
  df_populations = simulation$get_count_history() %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(Identity = paste0(genotype, epistate)) %>% 
    dplyr::rename(Generation = time, Population = count) %>% 
    dplyr::select(Generation, Identity, Population)
  
  # Tumour edges
  df_edges = simulation$get_lineage_graph() %>% 
    # dplyr::mutate(lblr = paste(ancestor, progeny), lbrl = paste(progeny, ancestor)) %>% 
    dplyr::distinct(ancestor, progeny) %>% 
    dplyr::rename(Parent = ancestor, Identity = progeny) %>% 
    dplyr::select(Parent, Identity)
  
  # Remove loops with a visit
  C = "Wild-type"
  Q = X = NULL
  repeat{
   Cx = C[1] 
   Q = c(Q, Cx)
   
   w_edg = df_edges %>% dplyr::filter(Parent == Cx, !(Identity %in% Q))
   
   X = dplyr::bind_rows(X, w_edg)
   
   C = C[-1]
   C = c(C, w_edg$Identity)

   if(length(C) == 0) break
  }
  
  df_edges = X
  
  max_tumour_size = df_populations %>%  
    dplyr::group_by(Generation) %>% 
    dplyr::summarise(Population = sum(Population)) %>% 
    dplyr::pull(Population) %>% 
    max()

  max_tumour_size = max_tumour_size * 1.05

  # Wild-type dynamics
  WT_dynamics = df_populations %>% 
    dplyr::group_by(Generation) %>% 
    dplyr::summarise(Population = sum(Population)) %>% 
    dplyr::mutate(Identity = "Wild-type", Population = max_tumour_size - Population)

  T_WT_dynamics = dplyr::bind_rows(
    WT_dynamics,
    df_populations
  )

  Muller_df = ggmuller::get_Muller_df(df_edges, T_WT_dynamics)
  
  ggmuller::Muller_pop_plot(Muller_df, add_legend = TRUE) + 
    my_theme() +
    ggplot2::guides(fill = ggplot2::guide_legend('Species')) +
    ggplot2::scale_fill_manual(values =
                                 c(`Wild-type` = 'gainsboro', get_species_colors(simulation))
    )
}

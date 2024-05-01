
#'Plot a vaf spectrum
#'
#' @description Plot a VAF spectrum where clusters are colored by stick
#'   membership or mutation SBS cause
#'
#' @param muts Mutation table where causes and/or stick membership of mutations
#'   are annotated.
#' @param vaf_cut Lower bound on VAF (defalt: 0.2).
#' @param color_by Specify whether to color mutations by stick membership or
#'   SBS cause. It can be "stick" or "SBS" (default: "stick").
#'
#' @return A `ggplot` plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation)
#' sim$update_tissue("Liver", 2e3, 2e3)
#'
#' sim$add_mutant(name = "A",
#'                growth_rates = 1,
#'                death_rates = 0)
#' sim$place_cell("A", 1e3, 1e3)
#' sim$run_up_to_size("A",1e4)
#'
#' sim$add_mutant(name = "B",
#'                growth_rates = 3.5,
#'                death_rates = 0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_size("B",1e4)
#'
#' bbox <- sim$search_sample(c("A" = 100,"B" = 100), 50, 50)
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#' forest = sim$get_samples_forest()
#'
#' m_engine = build_mutation_engine(setup_code = "demo")
#' m_engine$add_mutant(mutant_name = "A",
#'                     passenger_rates = c(SNV = 5e-8))
#' m_engine$add_mutant(mutant_name = "B",
#'                     passenger_rates = c(SNV = 5e-8))
#' m_engine$add_exposure(time = 0, c(SBS1 = 0.2,SBS5 = 0.8))
#'
#' phylo_forest = m_engine$place_mutations(forest, 100)
#'
#' labels = get_events_table(forest)
#' seq_results = simulate_seq(phylo_forest, coverage = 100,
#'                            write_SAM = F)
#' muts = muts_to_sticks(phylo_forest, seq_results, labels)
#'
#' plot_vaf(muts)

plot_vaf <- function(muts, vaf_cut = 0.02, color_by = "stick") {

  if (color_by == "stick") {
    if (! "obs" %in% colnames(muts)) {
      stop("Stick annotation missing")
    }
    cls <- ggsci::pal_simpsons()(muts %>%
                                   filter(.data$obs != "pre-neoplastic") %>%
                                   dplyr::pull(.data$obs) %>% unique() %>%
                                   length())
    names(cls) <- muts %>% filter(.data$obs != "pre-neoplastic")  %>%
      dplyr::arrange(dplyr::desc(.data$cell_id)) %>% dplyr::pull(obs) %>%
      unique()

    cls["Subclonal"] <- "gainsboro"
    cls["pre-neoplastic"] <- "salmon"
  } else if (color_by == "SBS") {
    cls <-  rRACES:::get_signatures_colors()
    muts <- muts %>% dplyr::mutate(obs = .data$causes)
  } else {
    stop("color_by should be 'stick' or 'SBS'")
  }

  muts <- muts %>% dplyr::select(-starts_with("normal"))

  dim <- grepl(x = colnames(muts), pattern = "VAF") %>% sum()

  ggplot_theme <- ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(.3, "cm"),
      panel.background = ggplot2::element_rect(fill = 'white')
    )

  if (dim == 1) {
    label <- colnames(muts)[grepl(x = colnames(muts), pattern = "VAF")]
    colnames(muts)[grepl(x = colnames(muts), pattern = "VAF")] <- "VAF"

    ggplot2::ggplot(muts %>% filter(.data$VAF > vaf_cut)) +
      ggplot2::geom_histogram(ggplot2::aes(x = .data$VAF, fill = .data$obs),
                              binwidth = 0.01) +
      ggplot2::scale_fill_manual(values = cls) + ggplot_theme +
      ggplot2::labs(title = "VAF distribution", x = label)
  } else {

    labels <- colnames(muts)[grepl(x = colnames(muts), pattern = "VAF")]
    colnames(muts)[grepl(x = colnames(muts),
                         pattern = "VAF")] <- c("VAF_1", "VAF_2")

    multivaf <- ggplot2::ggplot(muts %>% filter(.data$VAF_1 > vaf_cut |
                                                  .data$VAF_2 > vaf_cut)) +
      ggplot2::geom_point(ggplot2::aes(x = .data$VAF_1,
                                       y = .data$VAF_2, color = .data$obs), ) +
      ggplot2::scale_color_manual(values = cls) + ggplot_theme +
      ggplot2::labs(title = "VAF multivariate distribution",
                    x = labels[1], y = labels[2])

    vafn <- ggplot2::ggplot(muts %>% filter(.data$VAF_1 > vaf_cut)) +
      ggplot2::geom_histogram(ggplot2::aes(x = .data$VAF_1, fill = .data$obs),
                              binwidth = 0.01) + ggplot_theme +
      ggplot2::labs(title = paste0("marginal ", labels[1]), x = labels[1]) +
      ggplot2::scale_fill_manual(values = cls)

    vafp <- ggplot2::ggplot(muts %>% filter(.data$VAF_2 > vaf_cut)) +
      ggplot2::geom_histogram(ggplot2::aes(x = .data$VAF_2, fill = .data$obs),
                              binwidth = 0.01) + ggplot_theme +
      ggplot2::labs(title = paste0("marginal ", labels[2]), x = labels[2]) +
      ggplot2::scale_fill_manual(values = cls)

    list(multivaf, vafn, vafp)
  }
}

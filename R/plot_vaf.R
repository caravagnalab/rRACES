
#'Plot a vaf spectrum
#'
#' @description Plot a VAF spectrum where clusters are colored by stick
#'   membership or mutation SBS cause
#'
#' @param muts Mutation table where causes and/or stick membership of mutations
#'   are annotated.
#' @param vaf_cut Lower bound on VAF (default: 0.02).
#' @param colors A custom list of colors for any stick. If NULL a deafult palette is chosen.
#'
#' @return A `ggplot` plot.
#' @export
#'
#' @examples
#' sim <- new(Simulation)
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

plot_vaf = function(muts, colors = NULL, vaf_cut = 0.02){

  if(!is.null(muts$label)){

    if(!is.null(colors)){
      cls <- colors
    }else{
      cls_name <- muts %>% pull(label) %>% unique()
      cls <- RColorBrewer::brewer.pal(length(cls_name), "Dark2")
      names(cls) <- cls_name
      cls["Subclonal"] <- "gainsboro"
    }

  }

  muts <- muts %>% dplyr::select(-starts_with("normal"))

  dim <- grepl(x = colnames(muts),pattern = "VAF")  %>% sum()

  lab <- colnames(muts)[grepl(x = colnames(muts),pattern = "VAF")]

  marginals <- lapply(lab, function(l){

    y = muts
    colnames(y)[colnames(y) == l] = "VAF"

    plot <- ggplot2::ggplot(y %>% filter(VAF > vaf_cut)) +
      ggplot2::geom_histogram(ggplot2::aes(x = VAF), binwidth = 0.01) +
      my_theme() + ggplot2::labs(x = l)

    if(!is.null(muts$label)){

      plot <- plot +
        ggplot2::geom_histogram(ggplot2::aes(x = VAF, fill = label),
                                binwidth = 0.01) +
        ggplot2::scale_fill_manual(values = cls)

    }
    plot
  })

  if (dim > 1) {

    labs <- colnames(muts)[grepl(x = colnames(muts), pattern = "VAF")]

    g <- expand.grid(l1 = 1:length(labs),l2 = 1:length(labs)) %>%
        filter(l2 > l1) %>% rowwise() %>%
        mutate(lab1 = labs[l1], lab2 = labs[l2]) %>% ungroup()

    multi_plots <- lapply(1:nrow(g),function(i){
      y <- muts
      colnames(y)[colnames(y) %in% c(g$lab1[i],g$lab2[i])] <- c("VAF_1",
                                                                "VAF_2")

      multivaf <-  ggplot2::ggplot(y) +
        ggplot2::geom_point(ggplot2::aes(x = VAF_1,y = VAF_2)) +
        my_theme() + ggplot2::labs(x = g$lab1[i], y = g$lab2[i])

      if(!is.null(y$label)){

        multivaf <- multivaf +
          ggplot2::geom_point(ggplot2::aes(x = VAF_1, y = VAF_2,
                                           color = label)) +
            scale_color_manual(values = cls) +
          labs(x = g$lab1[i], y = g$lab2[i])

      }

      multivaf

    })

    plot <- list(marginals,multi_plots)
  } else {
    plot <- marginals
  }

  plot
}


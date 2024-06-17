#' Plot Marginals of Variant Allele Frequency (VAF)
#'
#' This function generates scatter plots showing the marginal distributions of
#' Variant Allele Frequency (VAF) for pairs of samples on a specific chromosome.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param chromosome A character specifying the chromosome to analyze (default:
#'   all the chromosomes in `seq_res`).
#' @param colour_by A character indicating whether to color the scatter points
#'   by "causes" or "classes" (default: "causes").
#' @param samples A character vector specifying the sample names to include
#'   in the plot. (default: NULL. It includes all samples except the
#'   "normal_sample").
#' @param cuts A numeric vector specifying the range of VAF values to include
#'   in the plot (default: `c(0, 1)`).
#' @return A list of ggplot2 objects showing scatter plots of VAF marginals
#'   for pairs of samples.
#' @seealso `plot_VAF_histogram()`, `plot_VAF()`
#' @export
#'
#' @examples
#' # set the seed of the random number generator
#' set.seed(0)
#'
#' sim <- SpatialSimulation()
#' sim$add_mutant(name = "A",
#'                growth_rates = 0.1,
#'                death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' # sampling tissue
#' n_w <- n_h <- 10
#' ncells <- 0.8 * n_w * n_h
#' bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
#' sim$sample_cells("SampleA", bbox$lower_corner, bbox$upper_corner)
#'
#' # adding second mutant
#' sim$add_mutant(name = "B",
#'                growth_rates = 0.3,
#'                death_rates = 0.0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_time(300)
#'
#' # sampling tissue again
#' bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
#' sim$sample_cells("SampleB", bbox$lower_corner, bbox$upper_corner)
#'
#' forest <- sim$get_samples_forest()
#'
#' # placing mutations
#' m_engine <- MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_mutant(mutant_name="B", passenger_rates=c(SNV=5e-9))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' # simulating sequencing
#' seq_results <- simulate_seq(phylo_forest, coverage = 10, write_SAM = F)
#'
#' library(dplyr)
#'
#' # filter germinal mutations and normal sample
#' f_seq <- seq_results %>% select(!starts_with("normal")) %>%
#'      filter(causes!="germinal")
#'
#' # plotting the VAR marginals
#' plot_VAF_marginals(f_seq)
#'
#' # plotting the VAF marginals and labelling it
#' plot_VAF_marginals(f_seq, labels = f_seq["causes"])
#'
#' # deleting the mutation engine directory
#' unlink('demo', recursive = T)
plot_VAF_marginals <- function(seq_res, chromosomes = NULL, samples = NULL,
                           labels = NULL, cuts = c(0, 1)) {

  data <- seq_to_long(seq_res)

  if (!is.null(labels)) {
    if (!is(labels, "data.frame")) {
        stop("The parameter \"labels\" must be a data frame when non-NULL.")
    }

    if (length(labels) != 1) {
        stop(paste0("The parameters \"labels\" must be a data frame ",
                    "with one column when non-NULL."))
    }

    if (nrow(labels) != nrow(seq_res)) {
        stop(paste0("The parameters \"seq_res\" and \"labels\"",
                    " must have the same number of rows."))
    }

    data["labels"] <- labels
    label_name <- names(labels)
  }

  chromosomes <- validate_chromosomes(seq_res, chromosomes)

  data <- data %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(sample_name != "normal_sample") %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  if (!is.null(samples)) {
    if (any(!samples %in% unique(data$sample_name))) {
      stop("Invalid sample name in samples parameter")
    }
    data <- data %>% dplyr::filter(sample_name %in% samples)
  }

  if (length(unique(data$sample_name))<2) {
    stop("At least two samples are required.")
  }

  combinations <- utils::combn(unique(data$sample_name), m = 2)

  lapply(1:ncol(combinations), function(i) {
    couple <- combinations[, i]
    d1 <- data %>% dplyr::filter(sample_name == couple[1]) #%>%
      #  dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))
    d2 <- data %>% dplyr::filter(sample_name == couple[2]) #%>%
      #  dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))

    djoin <- dplyr::full_join(d1, d2, by = c("chr", "from", "ref", "alt"))
    #djoin$col <- djoin[["labels.x"]]

    plot <- djoin %>% dplyr::filter(!is.na(sample_name.y))

    if (!is.null(labels)) {
        plot <- djoin %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = VAF.x, y = VAF.y,
                                                 col = labels.x))
    } else {
        plot <- djoin %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = VAF.x,
                                                 y = VAF.y))
    }

    plot +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::xlim(c(-0.01, 1.01)) +
      ggplot2::ylim(c(-0.01, 1.01)) +
      ggplot2::labs(x = couple[1], y = couple[2], col = "labels") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
  })
}

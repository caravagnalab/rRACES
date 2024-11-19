library(dplyr)

get_seq_data <- function(seq_res, sample, chromosomes)
{
  chromosomes <- validate_chromosomes(seq_res, chromosomes)

  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes)

  normal_data <- data %>% dplyr::filter(sample_name == "normal_sample")
  if (!sample %in% unique(data$sample_name)) {
    stop(paste("Inserted 'sample' is not available, available samples are:",
               paste0(unique(data$sample_name), collapse = ", ")))
  }
  tumour_data <- data %>% dplyr::filter(sample_name == sample)

  return(list(tumour=tumour_data, normal=normal_data))
}

#' Plot the genome-wide Depth Ratio (DR)
#'
#' @description This function plots the genome-wide Depth Ratio (DR)
#'   of a specific sample.
#'
#' @param seq_result A data frame containing sequencing results.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include
#'   in the plot (default: NULL, i.e., all the chromosomes).
#' @param N The number of mutations to sample for plotting (default: 5000).
#' @return A ggplot2 object showing the DR distribution across the genome.
#' @seealso `plot_VAF()`, `plot_BAF()`
#' @export
#'
#' @examples
#' # set the seed of the random number generator
#' set.seed(0)
#'
#' sim <- SpatialSimulation()
#' sim$add_mutant(name = "A",
#'                growth_rates = 0.2,
#'                death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' n_w <- n_h <- 10
#' ncells <- 0.8 * n_w * n_h
#' bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
#'
#' sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
#'
#' forest <- sim$get_samples_forest()
#' m_engine <- MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' seq_results <- simulate_seq(
#'   phylo_forest,
#'   chromosomes = c('22'),
#'   coverage = 10,
#'   write_SAM = F
#' )
#'
#' library(dplyr)
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>% filter(classes!="germinal")
#'
#' plot_DR(f_seq, sample="Sample")
#'
#' unlink('demo', recursive = T)
plot_DR <- function(
    seq_result,
    sample,
    chromosomes = NULL,
    N = 5000) {

  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_result) && ("mutations" %in% names(seq_result))) {

    # extract the field
    seq_res <- seq_result["mutations"]
  } else {
    seq_res <- seq_result
  }

  if (!("normal_sample.occurrences" %in% names(seq_res)
        && "normal_sample.coverage" %in% names(seq_res)
        && "normal_sample.VAF" %in% names(seq_res))) {
    stop(paste("The parameter \"seq_result\" does not contain",
               "the mandatory normal sample \"normal_sample\"."))
  }

  data <- get_seq_data(seq_res, sample, chromosomes)

  Ntotal <- nrow(data$tumour)
  N = min(N, Ntotal)

  d <- data$tumour %>%
    dplyr::left_join(data$normal, suffix = c(".tumour", ".normal"),
                     by = c("chr", "from", "ref", "alt")) %>%
    dplyr::mutate(DR = DP.tumour / DP.normal) %>%
    dplyr::sample_n(N) %>%
    dplyr::arrange(chr, from) %>%
    dplyr::mutate(abs_pos = 1:dplyr::n())

  chr_limits <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(abs_pos == min(abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))

  chr_means <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>%
    dplyr::pull(mean_pos)

  d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = abs_pos, y = DR)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid",
                        alpha = 0.4) +
    ggplot2::geom_point(size = 0.2, alpha = 0.3) +
    ggplot2::labs(x = NULL, y = "DR") +
    ggplot2::lims(y = c(0, NA)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means,
                                labels = unique(d$chr))
}


#' Plot the genome-wide B-Allele Frequency (BAF)
#'
#' @description This function plots the genome-wide B-Allele Frequency (BAF)
#'   of a specific sample.
#'
#' @param seq_result A data frame containing sequencing results.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include
#'   in the plot (default: NULL, i.e., all the chromosomes).
#' @param cuts A numeric vector specifying the range of BAF values to include
#'   in the plot (default: c(0, 1)).
#' @param N The number of mutations to sample for plotting (default: 5000).
#' @return A ggplot2 object showing the BAF distribution across the genome.
#' @seealso `plot_VAF()`, `plot_DR()`
#' @export
#'
#' @examples
#' # set the seed of the random number generator
#' set.seed(0)
#'
#' sim <- SpatialSimulation()
#' sim$add_mutant(name = "A",
#'                growth_rates = 0.2,
#'                death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' n_w <- n_h <- 10
#' ncells <- 0.8 * n_w * n_h
#' bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
#'
#' sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
#'
#' forest <- sim$get_samples_forest()
#' m_engine <- MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' seq_results <- simulate_seq(
#'   phylo_forest,
#'   chromosomes = c('22'),
#'   coverage = 10,
#'   write_SAM = F
#' )
#'
#' library(dplyr)
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>% filter(classes!="germinal")
#'
#' plot_BAF(f_seq, sample="Sample")
#'
#' unlink('demo', recursive = T)
plot_BAF <- function(
    seq_result,
    sample,
    chromosomes = NULL,
    cuts = c(0, 1),
    N = 5000) {

  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_result) && ("mutations" %in% names(seq_result))) {

    # extract the field
    seq_res <- seq_result["mutations"]
  } else {
    seq_res <- seq_result
  }

  data <- get_seq_data(seq_res, sample, chromosomes)

  Ntotal <- nrow(data$tumour)
  N = min(N, Ntotal)

  d <- data$tumour %>%
    #dplyr::group_by(chr) %>%
    dplyr::sample_n(N) %>%
    #dplyr::ungroup() %>%
    dplyr::arrange(chr, from) %>%
    dplyr::mutate(abs_pos = 1:dplyr::n()) %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  chr_limits <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(abs_pos == min(abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))

  chr_means <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>%
    dplyr::pull(mean_pos)

  d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = abs_pos, y = VAF)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid",
                        alpha = 0.4) +
    ggplot2::geom_point(size = 0.2, alpha = 0.3) +
    ggplot2::labs(x = NULL, y = "BAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr))
}


#' Plot the genome-wide Variant Allele Frequency (VAF)
#'
#' @description This function plots the genome-wide Variant Allele
#'   Frequency (VAF) of a specific sample.
#'
#' @param seq_result A data frame containing sequencing results.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include
#'   in the plot (default: NULL, i.e., all the chromosomes).
#' @param cuts A numeric vector specifying the range of VAF values to include
#'   in the plot (default: c(0, 1)).
#' @param N The number of mutations to sample for plotting (default: 5000).
#' @return A ggplot2 object showing the VAF distribution across the genome.
#' @seealso `plot_BAF()`, `plot_DR()`, `plot_VAF_histogram()`
#' @export
#'
#' @examples
#' # set the seed of the random number generator
#' set.seed(0)
#'
#' sim <- SpatialSimulation()
#' sim$add_mutant(name = "A",
#'                growth_rates = 0.2,
#'                death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' n_w <- n_h <- 10
#' ncells <- 0.8 * n_w * n_h
#' bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
#'
#' sim$sample_cells("Sample", bbox$lower_corner, bbox$upper_corner)
#'
#' forest <- sim$get_samples_forest()
#' m_engine <- MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' seq_results <- simulate_seq(
#'   phylo_forest,
#'   chromosomes = c('22'),
#'   coverage = 10,
#'   write_SAM = F
#' )
#'
#' library(dplyr)
#'
#' # filter germinal mutations
#' f_seq <- seq_results$mutations %>% filter(classes!="germinal")
#'
#' plot_VAF(f_seq, sample="Sample")
#'
#' unlink('demo', recursive = T)
plot_VAF <- function(seq_result, sample,
    chromosomes = NULL,
    cuts = c(0, 1), N = 5000) {

  # if the type of seq_res is a list and seq_res contains a field "mutations"
  if (is.list(seq_result) && ("mutations" %in% names(seq_result))) {

    # extract the field
    seq_res <- seq_result["mutations"]
  } else {
    seq_res <- seq_result
  }

  data <- get_seq_data(seq_res, sample, chromosomes)

  Ntotal = nrow(data$tumour)
  N = min(N, Ntotal)

  d <- data$tumour %>%
    dplyr::arrange(chr, from) %>%
    #dplyr::group_by(chr) %>%
    dplyr::sample_n(N) %>%
    #dplyr::ungroup() %>%
    dplyr::mutate(abs_pos = 1:dplyr::n()) %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))

  chr_limits <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(abs_pos == min(abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))

  chr_means <- d %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>%
    dplyr::pull(mean_pos)

  d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = abs_pos, y = VAF)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid",
                        alpha = 0.4) +
    ggplot2::geom_point(size = 0.5, alpha = 0.4) +
    ggplot2::labs(x = NULL, y = "VAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr))
}

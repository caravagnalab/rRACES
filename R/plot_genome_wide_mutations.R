#' Plot Depth Ratio (DR) Genome-wide
#'
#' This function generates a plot showing the Depth Ratio (DR) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot. Default is all autosomes and sex chromosomes.
#' @param N The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the DR distribution across the genome.
#' @export
#' 
#' @examples
#' sim <- new(Simulation)
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
#' m_engine <- build_mutation_engine(setup_code = "demo")
#' 
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#' 
#' phylo_forest <- m_engine$place_mutations(forest, 10)
#' 
#' seq_results <- simulate_seq(
#'   phylo_forest, 
#'   chromosomes = c('22'), 
#'   coverage = 10, 
#'   write_SAM = F
#' )
#' 
#' plot_DR_gw(seq_results, sample="Sample", chromosomes=c("22"))
#' 
#' unlink('demo', recursive = T)
plot_DR_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), N = 5000) {
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosome passed as input")
  
  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes == "germinal") %>%
    dplyr::mutate(mut_id = paste(chr, from, to, sep = ":"))
  
  normal_data <- data %>% dplyr::filter(sample_name == "normal_sample")
  if (!sample %in% unique(data$sample_name)) {
    stop(paste("Inserted 'sample' is not available, available samples are:",
               paste0(unique(data$sample_name), collapse = ", ")))
  }
  tumour_data <- data %>% dplyr::filter(sample_name == sample) #%>% dplyr::filter(VAF <= max(cuts), VAF >= min(cuts)) 
  
  Ntotal <- nrow(tumour_data)
  N = min(N, Ntotal)
  
  d <- tumour_data %>%
    dplyr::left_join(normal_data, suffix = c(".tumour", ".normal"), by = "mut_id") %>%
    dplyr::mutate(DR = DP.tumour / DP.normal) %>%
    #dplyr::group_by(chr.tumour) %>%
    dplyr::sample_n(N) %>%
    #dplyr::ungroup() %>%
    dplyr::arrange(chr.tumour, from.tumour) %>%
    dplyr::mutate(abs_pos = 1:dplyr::n())
  
  chr_limits <- d %>%
    dplyr::group_by(chr.tumour) %>%
    dplyr::filter(abs_pos == min(abs_pos)) %>%
    dplyr::pull(abs_pos)
  chr_limits <- c(chr_limits, max(d$abs_pos))
  
  chr_means <- d %>%
    dplyr::group_by(chr.tumour) %>%
    dplyr::summarise(mean_pos = (max(abs_pos) + min(abs_pos)) / 2) %>%
    dplyr::pull(mean_pos)
  
  d %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = abs_pos, y = DR)) +
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed", alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid", alpha = 0.4) +
    ggplot2::geom_point(size = 0.2, alpha = 0.3) +
    ggplot2::labs(x = "", y = "DR") +
    ggplot2::lims(y = c(0, NA)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr.tumour)) +
    #ggplot2::geom_hline(yintercept = stats::median(d$DR), col = "indianred", linetype = "dashed", linewidth = 1) +
    ggplot2::ggtitle(sample, subtitle = paste0('N = ', N, ' (', round(100 * N / Ntotal),'%)'))
}


#' Plot B-Allele Frequency (BAF) Genome-wide
#'
#' This function generates a plot showing the B-Allele Frequency (BAF) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot.
#'                    Default is all autosomes and sex chromosomes.
#' @param cuts A numeric vector specifying the range of BAF values to include in the plot.
#'             Default is c(0, 1).
#' @param N The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the BAF distribution across the genome.
#' @export
#' 
#' @examples
#' sim <- new(Simulation)
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
#' m_engine <- build_mutation_engine(setup_code = "demo")
#' 
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#' 
#' phylo_forest <- m_engine$place_mutations(forest, 10)
#' 
#' seq_results <- simulate_seq(
#'   phylo_forest, 
#'   chromosomes = c('22'), 
#'   coverage = 10, 
#'   write_SAM = F
#' )
#' 
#' plot_BAF_gw(seq_results, sample="Sample", chromosomes=c("22"))
#' 
#' unlink('demo', recursive = T)
plot_BAF_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), cuts = c(0, 1), N = 5000) {
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosome passed as input")
  
  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes == "germinal")
  
  if (!sample %in% unique(data$sample_name)) {
    stop(paste("Inserted 'sample' is not available, available samples are:",
               paste0(unique(data$sample_name), collapse = ", ")))
  }
  tumour_data <- data %>% dplyr::filter(sample_name == sample)
  
  Ntotal <- nrow(tumour_data)
  N = min(N, Ntotal)
  
  d <- tumour_data %>%
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
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed", alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid", alpha = 0.4) +
    ggplot2::geom_point(size = 0.2, alpha = 0.3) +
    ggplot2::labs(x = "", y = "BAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr)) +
    #ggplot2::geom_hline(yintercept = stats::median(d$VAF), col = "indianred", linetype = "dashed", linewidth = 1) +
    ggplot2::ggtitle(sample, subtitle = paste0('N = ', N, ' (', round(100 * N / Ntotal),'%)'))
}


#' Plot Variant Allele Frequency (VAF) Genome-wide
#'
#' This function generates a plot showing the Variant Allele Frequency (VAF) genome-wide for a specific sample.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param sample The name of the sample for which the plot will be generated.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot.
#'                    Default is all autosomes and sex chromosomes.
#' @param cuts A numeric vector specifying the range of VAF values to include in the plot.
#'             Default is c(0, 1).
#' @param N The number of mutations to sample for plotting. Default is 5000.
#' @return A ggplot2 object showing the VAF distribution across the genome.
#' @export
#' 
#' @examples
#' sim <- new(Simulation)
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
#' m_engine <- build_mutation_engine(setup_code = "demo")
#' 
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#' 
#' phylo_forest <- m_engine$place_mutations(forest, 10)
#' 
#' seq_results <- simulate_seq(
#'   phylo_forest, 
#'   chromosomes = c('22'), 
#'   coverage = 10, 
#'   write_SAM = F
#' )
#' 
#' plot_VAF_gw(seq_results, sample="Sample", chromosomes=c("22"))
#' 
#' unlink('demo', recursive = T)
plot_VAF_gw <- function(seq_res, sample, chromosomes = paste0(c(1:22, "X", "Y")), cuts = c(0, 1), N = 5000) {
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosome passed as input")
  
  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes != "germinal")
  
  if (!sample %in% unique(data$sample_name)) {
    stop(paste("Inserted 'sample' is not available, available samples are:",
               paste0(unique(data$sample_name), collapse = ", ")))
  }
  tumour_data <- data %>% dplyr::filter(sample_name == sample)
  
  Ntotal = nrow(tumour_data)
  N = min(N, Ntotal)
  
  d <- tumour_data %>%
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
    ggplot2::geom_vline(xintercept = chr_means, linetype = "dashed", alpha = 0.2) +
    ggplot2::geom_vline(xintercept = chr_limits, linetype = "solid", alpha = 0.4) +
    ggplot2::geom_point(size = 0.5, alpha = 0.4) +
    ggplot2::labs(x = "", y = "VAF") +
    ggplot2::lims(y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = chr_means, labels = unique(d$chr)) +
    #ggplot2::geom_hline(yintercept = stats::median(d$VAF), col = "indianred", linetype = "dashed", linewidth = 1) +
    ggplot2::ggtitle(sample) +
    ggplot2::ggtitle(sample, subtitle = paste0('N = ', N, ' (', round(100 * N / Ntotal),'%)'))
}

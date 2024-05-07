#' Plot Histogram of Variant Allele Frequency (VAF)
#'
#' This function generates a histogram showing the distribution of Variant Allele Frequency (VAF) across samples and chromosomes.
#'
#' @param seq_res A data frame containing sequencing results in wide format.
#' @param chromosomes A character vector specifying the chromosomes to include in the plot.
#'                    Default is all autosomes and sex chromosomes.
#' @param samples A character vector specifying the sample names to include in the plot.
#'                Default is NULL, which includes all samples except the "normal_sample".
#' @param colour_by A character indicating whether to color the histogram bars by "causes" or "classes".
#'                  Default is "causes".
#' @param cuts A numeric vector specifying the range of VAF values to include in the plot.
#'             Default is c(0, 1).
#' @return A ggplot2 object showing the histogram of VAF.
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
#' # sampling tissue 
#' n_w <- n_h <- 10
#' ncells <- 0.8 * n_w * n_h
#' bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
#' sim$sample_cells("SampleA", bbox$lower_corner, bbox$upper_corner)
#' 
#' # adding second mutant
#' sim$add_mutant(name = "B",
#'                growth_rates = 0.5,
#'                death_rates = 0.0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_time(200)   
#' 
#' # sampling tissue again
#' bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
#' sim$sample_cells("SampleB", bbox$lower_corner, bbox$upper_corner)
#' 
#' forest <- sim$get_samples_forest()
#' 
#' # placing mutations
#' m_engine <- build_mutation_engine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8))
#' m_engine$add_mutant(mutant_name="B", passenger_rates=c(SNV=5e-9))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#' 
#' phylo_forest <- m_engine$place_mutations(forest, 10)
#' 
#' # simulating sequencing 
#' seq_results <- simulate_seq(
#'   phylo_forest, 
#'   chromosomes = c('22'), 
#'   coverage = 10, 
#'   write_SAM = F
#' )
#' 
#' # plotting histogram of the VAF 
#' plot_histogram_vaf(seq_results, colour_by="causes", chromosomes=c("22"))
#' 
#' # deleting the mutation engine directory
#' unlink('demo', recursive = T)
plot_histogram_vaf <- function(
    seq_res,
    chromosomes = paste0(c(1:22, "X", "Y")),
    samples = NULL,
    colour_by = "causes",
    cuts = c(0, 1)
) {
  if (!(colour_by %in% c("classes", "causes")))
    stop("Colour_by parameter can be either 'causes' or 'classes'")
  if (any(!chromosomes %in% paste0(c(1:22, "X", "Y"))))
    stop("Invalid chromosomes passed as input")
  
  data <- seq_to_long(seq_res) %>%
    dplyr::mutate(chr = factor(chr, levels = chromosomes)) %>%
    dplyr::filter(chr %in% chromosomes) %>%
    dplyr::filter(classes != "germinal") %>%
    dplyr::filter(sample_name != "normal_sample") %>%
    dplyr::filter(VAF <= max(cuts), VAF >= min(cuts))
  
  if (!is.null(samples)) {
    if (any(!samples %in% unique(data$sample_name)))
      stop("Invalid sample name in samples parameter")
    data <- data %>% dplyr::filter(sample_name %in% samples)
  }
  
  data$col <- data[[colour_by]]
  
  data %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = VAF, fill = col)) +
    ggplot2::geom_histogram(binwidth=0.01, alpha=0.5, position = 'identity') +
    ggplot2::xlim(x = c(-0.01, 1.01)) +
    ggplot2::facet_grid(sample_name ~ chr, scales = 'free_y') +
    ggplot2::theme_bw() +
    ggplot2::labs(col = colour_by, fill = colour_by) +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    ggplot2::theme(legend.position = "bottom") #axis.text.x = element_text(size = 5)
}

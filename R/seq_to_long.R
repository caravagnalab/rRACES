#' Convert Sequencing Results from Wide to Long Format
#'
#' This function takes sequencing results in wide format and converts them into
#' a long format data frame. It extracts sample names from column names,
#' processes each sample separately, and then binds them together. Finally, it
#' renames and reorders columns to match the desired output format.
#'
#' @param seq_results A data frame containing sequencing results in wide format.
#' @return A data frame in long format with columns "`chr`", "`from`", "`to`",
#' "`ref`", "`alt`", "`NV`", "`DP`", "`VAF`", and "`sample_name`".
#' @export
#'
#' @examples
#' # Example data frame in wide format
#' seq_results <- data.frame(chr = c("chr1", "chr2"),
#'                           chr_pos = c(100, 200),
#'                           ref = c("A", "C"),
#'                           alt = c("T", "G"),
#'                           causes = c("SBS5", "SBS1"),
#'                           classes = c("germinal", "passneger"),
#'                           Sample.A.occurrences = c(10, 90),
#'                           Sample.A.coverage = c(100, 100),
#'                           Sample.A.VAF = c(0.1, 0.9),
#'                           normal_sample.occurrences = c(45, 52),
#'                           normal_sample.coverage = c(100, 100),
#'                           normal_sample.VAF = c(0.45, 0.52))
#' seq_results
#'
#' # Convert to long format
#' seq_to_long(seq_results)
seq_to_long <- function(seq_results) {
  # Extract sample names from column names
  sample_names <- strsplit(colnames(seq_results)[grepl(".VAF",
                                                       colnames(seq_results),
                                                       fixed = TRUE)],
                           ".VAF") %>% unlist()

  # Process each sample separately to create a list of data frames
  seq_df <- lapply(sample_names, function(sn) {
    # Select relevant columns for the current sample
    cc <- c("chr", "chr_pos", "ref", "alt", "causes", "classes",
            colnames(seq_results)[grepl(paste0(sn, "."),
                                        colnames(seq_results), fixed = TRUE)])

    # Rename columns and add sample_name column
    seq_results[, cc] %>%
      `colnames<-`(c("chr", "chr_pos", "ref", "alt", "causes", "classes",
                     "occurences", "coverage", "VAF")) %>%
      dplyr::mutate(sample_name = sn)
  }) %>% do.call("bind_rows", .)

  # Rename and reorder columns
  seq_df %>%
    dplyr::rename(from = chr_pos, DP = coverage,
                  NV = occurences) %>%
    dplyr::mutate(to = from)
}
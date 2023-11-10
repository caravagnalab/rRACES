#' @importFrom rlang .data
setMethod("show", "Rcpp_Simulation", function(object) {
  sim_status <- (length(object$get_species())
                 & nrow(object$get_cells()))
  sam_status <- FALSE
  mut_status <- FALSE

  sim_status <- ifelse(
    sim_status,
    crayon::bgGreen(crayon::white(" D ")),
    crayon::bgRed(crayon::white(" D "))
  )

  sam_status <- ifelse(
    sam_status,
    crayon::bgGreen(crayon::white(" S ")),
    crayon::bgRed(crayon::white(" S "))
  )

  mut_status <- ifelse(
    mut_status,
    crayon::bgGreen(crayon::white(" M ")),
    crayon::bgRed(crayon::white(" M "))
  )

  cli::cli_rule(
    left = paste(
      crayon::bgYellow(" rRACES "),
      sim_status, sam_status, mut_status,
      object$get_name()
    ),
    right = paste0(
      crayon::red("\u25A3 "),
      object$get_tissue_name(),
      " [",
      paste(object$get_tissue_size(), collapse = "x"),
      "]",
      crayon::red("  \u23f1 "),
      signif(object$get_clock(), digits = 3)
    )
  )

  species <- object$get_species()
  if (length(species) == 0) {
    cat('\n')
    cli::cli_alert_danger("The simulation has no species yet!")
  } else {
    counts_tab <- object$get_counts()

    my_tab <- species %>%
      dplyr::left_join(counts_tab, by = c("genotype", "epistate")) %>%
      dplyr::mutate(species = paste0(.data$genotype, .data$epistate)) %>%
      #dplyr::select(.data$species, ancestor, time, .data$growth_rate,
      dplyr::select(.data$species, .data$growth_rate,
                    .data$death_rate, .data$switch_rate, .data$counts)
    colnames(my_tab)[2:4] <- c(" \u03BB ", " \u03B4 ", " \u03B5 ")

    cat("   ",
        knitr::kable(my_tab, format = "rst", align = "rcrccc"),
        sep = "\n   ")
  }
})
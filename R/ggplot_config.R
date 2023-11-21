my_theme <- function() {
  ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

get_colors_for <- function(df, column_name, pal_name) {
  # Unpaired
  colors <- NULL
  if (nrow(df) > 0) {
    num_of_values <- df[, column_name] %>% length()

    if (num_of_values < 3) {
      colors <- RColorBrewer::brewer.pal(3, pal_name)

      colors <- colors[seq_len(num_of_values)]
    } else {
      colors <- RColorBrewer::brewer.pal(num_of_values, pal_name)
    }

    names(colors) <- df[, column_name]
  }

  return(colors)
}

get_species_colors <- function(species) {

  paired_species <- species %>%
    dplyr::filter(!is.na(.data$epistate)) %>%
    dplyr::mutate(
      species = paste0(.data$genotype, .data$epistate)
    ) %>%
    dplyr::arrange(.data$genotype)

  unpaired_species <- species %>%
    dplyr::filter(is.na(.data$epistate)) %>%
    dplyr::mutate(
      species = paste0(.data$genotype, .data$epistate)
    ) %>%
    dplyr::arrange(.data$genotype)

  # Unpaired
  unp_colour <- get_colors_for(unpaired_species, "genotype", "Dark2")

  # Paired
  pai_colour <- get_colors_for(paired_species, "species", "Paired")

  return(c(pai_colour, unp_colour))
}





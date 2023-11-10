my_theme <- function() {
  ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

get_species_colors <- function(x) {
  
  paired_species = x$get_species() %>% 
    dplyr::filter(!is.na(.data$epistate)) %>% 
    dplyr::mutate(
      species = paste0(.data$genotype, .data$epistate)
    ) %>% 
    dplyr::arrange(.data$genotype)
  
  unpaired_species = x$get_species() %>% 
    dplyr::filter(is.na(.data$epistate)) %>% 
    dplyr::mutate(
      species = paste0(.data$genotype, .data$epistate)
    ) %>% 
    dplyr::arrange(.data$genotype)
  
  # Unpaired
  unp_colour = NULL
  if(nrow(unpaired_species) > 0){
    unp_colour = RColorBrewer::brewer.pal(unpaired_species$genotype %>% length(), "Dark2")
    
    if (nrow(unpaired_species) < length(unp_colour)) {
      unp_colour <- unp_colour[1:nrow(unpaired_species)]
    }
    
    names(unp_colour) <- unpaired_species$genotype
  }
  
  # Paired
  pai_colour = NULL
  if(nrow(paired_species) > 0){
    pai_colour = RColorBrewer::brewer.pal(paired_species$species %>% length(), "Paired")
    
    if (nrow(paired_species) < length(pai_colour)) {
      pai_colour <- pai_colour[1:nrow(paired_species)]
    }
    
    names(pai_colour) <- paired_species$species
  }

  return(c(pai_colour, unp_colour))
}





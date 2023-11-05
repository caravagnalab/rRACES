# My own private getters for the species table
ExistsSpecies = function(x, genotype){
  # x$species %>% 
  #   dplyr::filter(genotype == !!genotype) %>% 
  #   nrow()
}

AncestorOf = function(x, genotype){
  x$species %>% 
    dplyr::filter(genotype == !!genotype) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::pull(.data$ancestor)
}

TimeOf = function(x, genotype){
  x$species %>% 
    dplyr::filter(genotype == !!genotype) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::pull(.data$time)
}

GRateOf = function(x, genotype, epistate) {
  x$species %>%
    dplyr::filter(genotype == !!genotype, epistate == !!epistate) %>%
    dplyr::pull(.data$rgrowth)
}

DRateOf = function(x, genotype, epistate) {
  x$species %>%
    dplyr::filter(genotype == !!genotype, epistate == !!epistate) %>%
    dplyr::pull(.data$rdeath)
}

ERateOf = function(x, genotype, epistate) {
  x$species %>%
    dplyr::filter(genotype == !!genotype, epistate == !!epistate) %>%
    dplyr::pull(.data$repigenetic)
}
# 
# obj_status = function(x){
#   can_simulate
#   if(x$has_species & x$has_initial_cell)
#   "rRACES object",
#   "Species added",
#   "Initial state"
#   
# }
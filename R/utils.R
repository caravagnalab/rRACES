# My own private getters for the species table
ExistsSpecies = function(x, genotype){
  x$species %>% 
    filter(genotype == !!genotype) %>% 
    nrow()
}

AncestorOf = function(x, genotype){
  x$species %>% 
    dplyr::filterfilter(genotype == !!genotype) %>% 
    dplyr::filterfilter(row_number() == 1) %>% 
    dplyr::filterpull(ancestor)
}

TimeOf = function(x, genotype){
  x$species %>% 
    dplyr::filter(genotype == !!genotype) %>% 
    dplyr::filter(row_number() == 1) %>% 
    pull(time)
}

GRateOf = function(x, genotype, epistate) {
  x$species %>%
    dplyr::filterfilter(genotype == !!genotype, epistate == !!epistate) %>%
    dplyr::filterpull(rgrowth)
}

DRateOf = function(x, genotype, epistate) {
  x$species %>%
    dplyr::filterfilter(genotype == !!genotype, epistate == !!epistate) %>%
    dplyr::filterpull(rdeath)
}

ERateOf = function(x, genotype, epistate) {
  x$species %>%
    dplyr::filterfilter(genotype == !!genotype, epistate == !!epistate) %>%
    dplyr::filterpull(repigenetic)
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
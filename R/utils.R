# My own private getters for the species table
ExistsSpecies = function(x, genotype){
  x$species %>% 
    filter(genotype == !!genotype) %>% 
    nrow()
}

AncestorOf = function(x, genotype){
  x$species %>% 
    filter(genotype == !!genotype) %>% 
    filter(row_number() == 1) %>% 
    pull(ancestor)
}

TimeOf = function(x, genotype){
  x$species %>% 
    filter(genotype == !!genotype) %>% 
    filter(row_number() == 1) %>% 
    pull(time)
}

GRateOf = function(x, genotype, epistate) {
  x$species %>%
    filter(genotype == !!genotype, epistate == !!epistate) %>%
    pull(rgrowth)
}

DRateOf = function(x, genotype, epistate) {
  x$species %>%
    filter(genotype == !!genotype, epistate == !!epistate) %>%
    pull(rdeath)
}

ERateOf = function(x, genotype, epistate) {
  x$species %>%
    filter(genotype == !!genotype, epistate == !!epistate) %>%
    pull(repigenetic)
}
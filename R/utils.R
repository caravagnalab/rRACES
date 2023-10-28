# My own private getters for the species table
AncestorOf = function(x,s){
  x$species %>% 
    filter(species == s) %>% 
    filter(row_number() == 1) %>% 
    pull(ancestor)
}

TimeOf = function(x,s){
  x$species %>% 
    filter(species == s) %>% 
    filter(row_number() == 1) %>% 
    pull(time)
}

GRateOf = function(x, s, e) {
  x$species %>%
    filter(species == s, epistate == e) %>%
    pull(rgrowth)
}

DRateOf = function(x,s,e){
  x$species %>%
    filter(species == s, epistate == e) %>%
    pull(rdeath)
}

ERateOf = function(x,s,e){
  x$species %>%
    filter(species == s, epistate == e) %>%
    pull(repigenetic)
}
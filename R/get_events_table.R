
#' Get stick table
#'
#' @description
#' Get subset of forest nodes corresponding to phylogenetic branches containing relevant biological events (drivers)
#'
#' @param forest The original forest object
#' has been derived
#'
#' @return A `tidyverse` tibble object
#' @export
#'
#' @examples
#' sim <- new(Simulation)
#' sim$add_mutant(name = "A", growth_rates = 0.08, death_rates = 0.01)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(60)
#' sim$sample_cells("MySample", c(500, 500), c(510, 510))
#' forest = sim$get_samples_forest()
#' forest$get_samples_info()
#' get_events_table(forest)


get_events_table = function(forest){
  
  sticks = forest$get_sticks()
  
  events_table = lapply(1:length(sticks), function(i){
    
    gen = forest$get_nodes() %>% as_tibble() %>% 
      filter(cell_id == sticks[[i]][length(sticks[[i]])]) %>% pull(mutant)
    
    epi = forest$get_nodes() %>% as_tibble() %>% 
      filter(cell_id == sticks[[i]][length(sticks[[i]])]) %>% pull(epistate)
   
    if (epi != "") {
      change = paste0(ifelse(epi == "+", "-", "+"), " -> ", epi)
      tibble(cell_id = sort(sticks[[i]][2:length(sticks[[i]])]),
             obs = paste0(gen, " ", change))
    } else{
      tibble(cell_id = sort(sticks[[i]][2:length(sticks[[i]])]), obs = gen
      )
    }
    
 }) %>% bind_rows()
  
  as_tibble(forest$get_nodes()) %>% full_join(events_table, by = "cell_id") %>% 
    mutate(obs = ifelse(is.na(obs),"Subclonal",obs)) %>% unique()
  
}





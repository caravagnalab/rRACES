
#' Get mutations-nodes map
#'
#' @description
#' Map mutations to nodes of the phylogenetic forest where they occured first.
#'
#' @param phylo_forest The original forest object where mutations has been
#'   placed.
#' @param seq_results Table annotating mutations genome coordinates and causes
#'   (the output of `simulate_seq()`).
#' @param cell_labels Table annotating the sticks (it can be the output of `get_events_table()`).
#' @param mutation_labels Table annotating user defined labels of the mutations.
#'
#' @return A `tidyverse` tibble object
#' @export
#'
#' @examples
#' library(rRACES)
#' library(dplyr)
#' sim <- new(Simulation)
#' sim$update_tissue("Liver", 2e3, 2e3)
#'
#' sim$add_mutant(name = "A",
#'                growth_rates = 1,
#'                death_rates = 0)
#' sim$place_cell("A", 1e3, 1e3)
#' sim$run_up_to_size("A",1e4)
#' sim$add_mutant(name = "B",
#'                growth_rates = 3.5,
#'                death_rates = 0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#'
#' sim$run_up_to_size("B",1e4)
#'
#' bbox <- sim$search_sample(c("A" = 100,"B" = 100), 50, 50)
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#' forest = sim$get_samples_forest()
#' labels = get_events_table(forest)
#'
#' m_engine = build_mutation_engine(setup_code = "demo")
#' m_engine$add_mutant(mutant_name = "A",
#'                     passenger_rates = c(SNV = 5e-8))
#' m_engine$add_mutant(mutant_name = "B",
#'                     passenger_rates = c(SNV = 5e-8))
#'
#' m_engine$add_exposure(time = 0, c(SBS1 = 0.2,SBS5 = 0.8))
#'
#' phylo_forest = m_engine$place_mutations(forest, 100)
#'
#' labels = get_events_table(forest)
#'
#' seq_results = simulate_seq(phylo_forest, coverage = 100, write_SAM = F)
#'
#' muts_to_sticks(phylo_forest, seq_results, cell_labels = labels)

muts_to_sticks = function(phylo_forest,seq_results,mutation_labels = NULL,cell_labels = NULL){
  
  sampled_snvs = phylo_forest$get_sampled_cell_mutations() %>% as_tibble() %>% 
    dplyr::select(-cell_id) %>% unique() %>% 
    mutate(id = paste0(
      "chr",
      chr,
      ":",
      chr_pos,
      ":",
      ref,
      ":",
      alt
    ))
  
  
  seq_results = seq_results %>% 
    mutate(id = paste0(
      "chr",
      chr,
      ":",
      chr_pos,
      ":",
      ref,
      ":",
      alt
    )) %>% as_tibble()
  
  if("classes" %in% colnames(seq_results)){
    
    seq_results = seq_results %>% filter(classes != "germinal")
    
  }
  
  seq_results = seq_results %>% as_tibble() %>% rowwise() %>%
    mutate(cell_id = phylo_forest$get_first_occurrences(SNV(
      chr, chr_pos, alt, ref
    ))[[1]]) %>%
    ungroup()
  
  if(!is.null(cell_labels)){
    map = full_join(seq_results,cell_labels, by = "cell_id") %>% 
      filter(!is.na(chr)) 
  }else if(!is.null(mutation_labels)){
    
    map = full_join(seq_results,mutation_labels, by = "id") %>% 
      filter(!is.na(chr)) 
    
  }else{
    
    map = seq_results
    
  }
  
  return(map)
  
}



#' Get mutations-nodes map
#'
#' @description
#' Map mutations to nodes of the phylogenetic forest where they occured first.
#'
#' @param phylo_forest The original forest object where mutations has been placed.
#' @param seq_results Table annotating mutations genome coordinates and causes (the output of simulate_seq).
#' @param labels Table annotating the sticks (output of get_events_table).
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
#' m_engine = build_mutation_engine(setup_code = "demo")
#' m_engine$add_mutant(mutant_name = "A",
#'                     passenger_rates = c(SNV = 5e-8))
#' m_engine$add_exposure(time = 0, c(SBS1 = 0.2,SBS5 = 0.8))
#' phylo_forest = m_engine$place_mutations(frel, 100)
#' labels = get_events_table(forest)
#' seq_results = simulate_seq(phylo_forest, coverage = 100,
#'                            epi_FACS = False, write_SAM = False)
#' muts_to_sticks(phylo_forest,seq_results,labels)

muts_to_sticks <- function(phylo_forest,seq_results,labels){

  sampled_mutations <- phylo_forest$get_sampled_cell_mutations() %>%
    as_tibble() %>%
    dplyr::select(-cell_id) %>% unique() %>%
    mutate(id = paste0("chr", chromosome, ":", chr_pos,
                       ":", ref, ":", alt))

  seq_results <- seq_results %>% 
    mutate(id = paste0("chr", chromosome, ":", chr_pos,
                       ":", ref, ":", alt)) %>%
    filter(classes != "germinal") %>% as_tibble()

  if(sum(grepl(x = colnames(seq_results), pattern = "P.VAF")) > 0){
    colnames(seq_results)[grep(x = colnames(seq_results), pattern = "P.VAF")] = "VAF_p"
    colnames(seq_results)[grep(x = colnames(seq_results), pattern = "N.VAF")] = "VAF_n"
    colnames(seq_results)[grep(x = colnames(seq_results), pattern = "P.coverage")] = "Depth_p" 
    colnames(seq_results)[grep(x = colnames(seq_results), pattern = "N.coverage")] = "Depth_n" 
    
  }
  
  
  seq_results = seq_results %>% as_tibble() %>% rowwise() %>%
    mutate(cell_id = phylo_forest$get_first_occurrences(SNV(
      chromosome, chr_pos, alt, ref
    ))[[1]]) %>%
    ungroup()
  
  
  map = full_join(seq_results,labels, by = "cell_id") %>% 
    filter(!is.na(chromosome)) %>% mutate(obs = 
                          ifelse(classes == "pre-neoplastic","pre-neoplastic",obs))
  
  return(map)
  
}

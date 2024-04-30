
#'Plot a vaf spectrum 
#'
#' @description
#' Plot a vaf spectrum where clusters are colored by stick membership or mutation SBS cause
#'
#' @param muts Mutation table where causes and/or stick membership of mutations are annotated.
#' @param vaf_cut Lower bound on vaf
#' @param color_by Specify whether to color mutations by stick membership or SBS cause. It can be "stick" or "SBS".
#'
#' @return A `ggplot` plot.
#' @export
#'
#' @examples
#'sim <- new(Simulation, "homogeneous_test",seed = sample(x = 1:5000,size = 1,replace = F),
#'            save_snapshot = F)
#' sim$duplicate_internal_cells <- T
#' sim$update_tissue("Liver", 2e3, 2e3)
#' sim$add_mutant(name = "A",
#'                growth_rates = 1,
#'                death_rates = 0)
#' sim$place_cell("A", 1e3, 1e3)
#' sim$run_up_to_size("A",1e4)
#' sim$add_mutant(name = "B",
#'                growth_rates = 3.5,
#'                death_rates = 0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_size("B",1e4)
#' bbox <- sim$search_sample(c("A" = 100,"B" = 100), 50, 50)
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#' forest = sim$get_samples_forest()
#' labels = get_events_table(forest)
#' m_engine = build_mutation_engine(setup_code = "demo")
#' m_engine$add_mutant(mutant_name = "A",
#' passenger_rates = c(SNV = 5e-8),
#' driver_SNVs = c(),
#' driver_CNAs = c())
#' m_engine$add_mutant(mutant_name = "B",
#' passenger_rates = c(SNV = 5e-8),
#' driver_SNVs = c(),
#' driver_CNAs = c())
#' m_engine$add_exposure(time = 0, c(SBS1 = 0.2,SBS5 = 0.8)) 
#' phylo_forest = m_engine$place_mutations(forest, 100)
#' labels = get_events_table(forest)
#' seq_results = simulate_seq(
#' phylo_forest,
#' coverage = 100,
#' epi_FACS =  F,
#' write_SAM = F)
#' muts = muts_to_sticks(phylo_forest,seq_results,labels)
#' plot_vaf(muts,vaf_cut = 0.02, color_by = "stick")


plot_vaf = function(muts,vaf_cut = 0.02, color_by = "stick"){

  if(color_by == "stick"){ 
    if(! "obs" %in% colnames(muts)){
      stop("Stick annotation missing")
    }
  cls = ggsci::pal_simpsons()(muts %>% 
                                filter(obs != "pre-neoplastic") %>% pull(obs) %>% unique() %>% length())
  names(cls) = muts %>% filter(obs != "pre-neoplastic")  %>% arrange(desc(cell_id)) %>% 
    pull(obs) %>% unique()
  cls["Subclonal"] = "gainsboro"
  cls["pre-neoplastic"] = "salmon"
  
  }else if(color_by == "SBS"){
    
  cls =  rRACES:::get_signatures_colors() 
  muts = muts %>% dplyr::mutate(obs = causes)
    
  }else{
    
     stop("color_by should be 'stick' or 'SBS'")
  }
  
  muts = muts %>% dplyr::select(-starts_with("normal"))
  
  dim = grepl(x = colnames(muts),pattern = "VAF")  %>% sum()  
  
  if(dim == 1){  
    
    label = colnames(muts)[grepl(x = colnames(muts),pattern = "VAF")]
    colnames(muts)[grepl(x = colnames(muts),pattern = "VAF")] = "VAF"
    
    ggplot(muts %>% filter(VAF > vaf_cut)) + 
      geom_histogram(aes(x = VAF, fill = obs), binwidth = 0.01) +
      scale_fill_manual(values = cls) + CNAqc:::my_ggplot_theme() + 
      labs(title = "VAF distribution", x = label)
    
  }else{
    
    labels = colnames(muts)[grepl(x = colnames(muts),pattern = "VAF")]
    colnames(muts)[grepl(x = colnames(muts),pattern = "VAF")] = c("VAF_1","VAF_2")
    
    multivaf =  ggplot(muts %>% filter(VAF_1 > vaf_cut | VAF_2 > vaf_cut)) + 
      geom_point(aes(x = VAF_1,y = VAF_2, color = obs),) +
      scale_color_manual(values = cls) + CNAqc:::my_ggplot_theme() + 
      labs(title = "VAF multivariate distribution",
           x = labels[1], y = labels[2])
    
    vafn = ggplot(muts %>% filter( VAF_1 > vaf_cut)) + 
      geom_histogram(aes(x = VAF_1, fill = obs),binwidth = 0.01) + CNAqc:::my_ggplot_theme() + 
      labs(title = paste0("marginal ",labels[1]),x = labels[1]) +
      scale_fill_manual(values = cls)
    
    vafp = ggplot(muts %>% filter( VAF_2 > vaf_cut)) + 
      geom_histogram(aes(x = VAF_2,fill = obs),binwidth = 0.01) + CNAqc:::my_ggplot_theme() + 
      labs(title = paste0("marginal ",labels[2]), x = labels[2]) +
      scale_fill_manual(values = cls)
    
    list(multivaf,vafn,vafp)
    
  }
  
  
}



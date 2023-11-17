#' Bounding box sampler
#' 
#' @description
#' From a simulation object, this function samples bounding boxes of 
#' a custom size until a number of cells of a certain species are not 
#' sampled. If this never happens, after a fixed number of attempts
#' the best box is returned.
#'
#' @param simulation A simulation object.
#' @param which The species name.
#' @param n The desired number of cells from species `which`.
#' @param n_w Width of the box.
#' @param n_h Height of the box
#'
#' @return coordinates for a bounding box.
#' @export
#'
#' @examples
#' sim <- new(Simulation, "bbox_sampler_test")
#' sim$add_genotype(name = "A", growth_rates = 0.08, death_rates = 0.01)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(60)
#' bbox = bbox_sampler(sim, "A", n = 15, n_w = 5, n_h = 5)
#' sim$sample_cells("A", bbox$p, bbox$q)
#' plot_tissue(sim)
bbox_sampler = function(simulation, which, n, n_w, n_h, nattempts = 100)
{
  # n_w = ceiling(sqrt(n))
  # n_h = ceiling(sqrt(n))
  
  t_size = simulation$get_tissue_size()
  
  sp <- cli::make_spinner()

  # where are cells of species which
  locations = simulation$get_cells() %>% 
    dplyr::mutate(species = paste0(.data$genotype, .data$epistate)) %>% 
    dplyr::filter(.data$species == which) %>% 
    dplyr::select(.data$position_x, .data$position_y)
  
  bof_p = bof_q = NULL
  bof_counts = -1
  
  repeat{
    sp$spin()
    
    location_id = sample(1:nrow(locations), 1)
    
    p = locations[location_id, ] %>% as.numeric()
    q = p + c(n_w, n_h)
    
    nc = simulation$get_cells(p, q) %>% 
      dplyr::mutate(species = paste0(.data$genotype, .data$epistate)) %>% 
      dplyr::filter(.data$species == which) %>% 
      nrow()
    
    if(nc > bof_counts)
    {
      bof_p = p
      bof_q = q
      bof_counts = nc
    }
    
    # cat(nc, '\n')
    
    if(nc >= n | nattempts == 0) break
    nattempts = nattempts - 1
  }
  
  sp$finish()
  
  return(list(p = p, q = q))
}

# bbox_sampler(sim, "A", 15)
run_simulation = function(
    genotypes,
    initial_cells,
    timed_events,
    tissue_name = "Custom Tissue",
    tissue_size = c(1200, 1200),
    death_activation_level = 150
    )
{
  command = '~/Downloads/RACES/drivers_sim'

  system(paste0(command, " ~/Downloads/RACES/examples/drivers_sim.json 200 -p"), intern = FALSE)
}


add_species = function(x,
                       name,
                       transition_rate = c(`+-` = 0.01, `-+` = 0.01),
                       growth_rate = c(`+` = 0.2, `-` = 0.15),
                       death_rate = c(`+` = 0.1, `-` = 0.1),
                       time_of_origin = 0,
                       ancestor = NULL,
                       driver = "KRAS G12D"
                       )
{
  if(is.null(x))
  {
    x = list()
    x$species = list()
    x$lineage = list()
  }

  nspecies = ifelse(
    is.null(x$species),
    0,
    x$species %>% length
    )

  if(nspecies == 0 & time_of_origin > 0)
  {
    cli::cli_abort("time_of_origin = 0 but no species, the first species needs to have 0")
  }

  if(nspecies == 0 & !is.null(ancestor))
  {
    cli::cli_abort("ancestor not NULL but no species, forcing WT to be the ancestor")
  }

  new_species = list(
    name = name,
    transition_rate = transition_rate,
    growth_rate = growth_rate,
    death_rate = death_rate
  )

  x$species = append(
    x$species,
    list(new_species)
  )

  new_lineage = list(
    parent = ancestor,
    child = name,
    time_of_origin = time_of_origin,
    driver = driver
  )

  x$lineage = append(
    x$lineage,
    list(new_lineage)
  )

  return(x)
}

x = add_species(NULL, name = "A") %>% add_species(name = "B")

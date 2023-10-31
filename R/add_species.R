#' Add a species to a simulation.
#' 
#' @description
#' This function adds one new species to the simulation. The species has its own
#' evolutionary parameter, an ancestor, and a time in which it is created. There
#' are special cases for the first species included, as well as several
#' consistency checks.
#' 
#' @param x A simulation.
#' @param name The genotype of the new species.
#' @param epigenetic_rates The epigenetic rates of the new species.
#' @param growth_rates The duplication rates of the new species.
#' @param death_rates The death rates of the new species.
#' @param ancestor The genotype from which the new species will arise.
#' @param time_of_origin The time at which the new species will arise.
#'
#' @return A simulation object.
#' @export
#'
#' @examples
#' create_simulation()
add_species = function(
    x,
    name = "A",
    epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
    growth_rates = c("+" = 0.2, "-" = 0.08),
    death_rates = c("+" = 0.1, "-" = 0.01),
    ancestor = NULL,
    time_of_origin = 0
)
{
  stopifnot(inherits(x, 'rraces'))
  
  # Check inputs -- first species
  if(!x$has_species){
    if(!is.null(time_of_origin) & time_of_origin != 0)
      cli::cli_abort("Time of origin has to be 0 for the first input species")
    if(!is.null(ancestor))
      cli::cli_abort("Ancestor has to be NULL for the first input species")
  }
  else
  {
    # Check inputs -- any other species
    if(name %in% x$species$species)
      cli::cli_abort("Species {.value {name}} already exists.")
    
    if(is.null(ancestor))
      cli::cli_abort("Ancestor cannot be NULL for this species")
    
    if(!(ancestor %in% x$species$species))
      cli::cli_abort("Ancestor {.value {ancestor}} does not exist.")
    
    t_ancestor = TimeOf(x, ancestor)
    
    if(t_ancestor >= time_of_origin)
      cli::cli_abort("Ancestor is born at time {.value {t_ancestor}}, so the species cannot be born at time {.field {time_of_origin}}")
  }
  
  # Check inputs -- common checks
  if(is.null(epigenetic_rates)) cli::cli_abort("Epigenetic rates are NULL")
  else{
    if(!all(names(epigenetic_rates) %in% c("+-", "-+"))) cli::cli_abort("Epigenetic rates are missing fields") 
    if(any(epigenetic_rates < 0)) cli::cli_abort("Epigenetic rates are negative") 
  }
  
  if(is.null(growth_rates)) cli::cli_abort("Growth rates are NULL")
  else{
    if(!all(names(growth_rates) %in% c("+", "-"))) cli::cli_abort("Growth rates are missing fields") 
    if(any(growth_rates < 0)) cli::cli_abort("Growth rates are negative") 
  }
  
  if(is.null(death_rates)) cli::cli_abort("Death rates are NULL")
  else{
    if(!all(names(death_rates) %in% c("+", "-"))) cli::cli_abort("Death rates are missing fields") 
    if(any(death_rates < 0)) cli::cli_abort("Death rates are negative") 
  }
  
  # We can add it now
  x$simulation$add_species(name = name,
                           epigenetic_rates = epigenetic_rates,
                           growth_rates = growth_rates,
                           death_rates = death_rates
  )
  
  x$species = x$species %>% 
    dplyr::bind_rows(
      dplyr::tibble(
        species = name,
        ancestor = ancestor,
        time = time_of_origin,
        epistate = '+',
        rgrowth = growth_rates['+'],
        rdeath = death_rates['+'],
        repigenetic = epigenetic_rates['+-']
      )
    ) %>% 
    dplyr::bind_rows(
      dplyr::tibble(
        species = name,
        ancestor = ancestor,
        time = time_of_origin,
        epistate = '-',
        rgrowth = growth_rates['-'],
        rdeath = death_rates['-'],
        repigenetic = epigenetic_rates['-+']
      )
    )
  
  x$has_species = TRUE
  
  cli::cli_alert_info("{crayon::underline('New species')} {.field {name}}, with ancestor {.field {ancestor}}")
  
  print(x)
  
  return(x)
}
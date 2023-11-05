#' Run a simulation.
#'
#' @description
#' This function runs a simulation forward in time, advancing the simulation to 
#' a specific point based on the provided parameters. The function takes a simulation
#' object `x` as the first argument and a list of named parameters (parameter `what`) 
#' that specify how to advance the simulation. 
#' 
#' The parameters can include:
#' - `time`: to move the simulation clock up to `time`;
#' - `genotype`/`epistate`/`size`: to move the simulation clock up to when the
#' species defined by `genotype` and `epistate` reaches size `size`;
#' - `genotype`/`epistate`/ `event`/`count`: to move the simulation clock up to 
#' when the model has fired `count` events of type `event` (any of `death`, `growth`
#' or `switch`) for the species defined by `genotype` and `epistate`.
#'  
#' The progress bar will show the percentage to completion of the appropriate task.
#'
#' @param x A simulation object.
#' @param what A list of named parameters that specify how to advance the simulation. 
#' The parameters are discussed in the description of this function.
#'
#' @return An updated simulation object.
#'
#' @export
#'
#' @examples
#'
#' # Example of using the function
#' x = create_simulation(output_dir="run_test")
#' x = add_genotype(x, "A")
#' x = set_initial_cell(x, "A", "+", c(50, 50))
#' x = run(x, list(time = 60))
#' x
run = function(
    x,
    what
    )
{
  # run_up_to_size
  # run_up_to_time
  # run_up_to_event
  
  stopifnot(inherits(x, 'rraces'))
  
  if(!x$has_species)
    cli::cli_abort("Species have to be added!")

  if(!x$has_initial_cell)
    cli::cli_abort("An initial cell has not been added!")
  
  # Time
  if('time' %in% names(what)) x = run_aux_time(x, what[['time']])
  else{
    if('size' %in% names(what)) x = run_aux_size(x, 
                                                 what[['genotype']], 
                                                 what[['epistate']], 
                                                 what[['size']])
    
    if('event' %in% names(what)) x = run_aux_event(x, 
                                                   event = what[['event']], 
                                                   genotype = what[['genotype']], 
                                                   epistate = what[['epistate']], 
                                                   count = what[['count']])
  }
  
  cli::cli_h3("New simulation state")
  
  print(x)
  
  x
}

run_aux_time = function(x, time)
{
  if(is.null(time))
    cli::cli_abort("NULL parameters!")
  
  if(!is.numeric(time))
    cli::cli_abort("Time requires a numeric.")
  
  if(time < 0)
    cli::cli_abort("Time cannot be negative.")
  
  if(x$simulation$get_clock() >= time)
    cli::cli_abort("Time smaller than the current clock")
  
  cli::cli_alert_info("Current clock {.field {x$simulation$get_clock()}}, target {.field {time}}.")
  
  x$simulation$run_up_to_time(time)
  
  return(x)
}

run_aux_size = function(x, genotype, epistate, size)
{
  if(is.null(genotype) | is.null(epistate) | is.null(size))
    cli::cli_abort("NULL parameters!")
  
  if(!is.numeric(size))
    cli::cli_abort("Size requires a numeric.")
  
  if(size < 0)
    cli::cli_abort("Size cannot be negative.")

  if(!(genotype %in% x$species$genotype))
    cli::cli_abort("Genotype does not exist in the simulation.")
  
  cc = x$simulation$get_counts() %>% 
    dplyr::filter(genotype == !!genotype, epistate == !!epistate) %>% 
    dplyr::pull(counts)
  
  cli::cli_alert_info("{.field {genotype}} with epistate {.field {epistate}} has count {.field {cc}}, target is {.field {size}}.")
  
  x$simulation$run_up_to_size(paste0(genotype, epistate), size)

  return(x)
}

run_aux_event = function(x, event, genotype, epistate, count)
{
  # sim$run_up_to_event("death", "B-", 100)
  
  if(is.null(event) | is.null(genotype) | is.null(epistate) | is.null(count))
    cli::cli_abort("NULL parameters!")
  
  if(!is.numeric(count))
    cli::cli_abort("Count requires a numeric.")
  
  if(count < 0)
    cli::cli_abort("Count cannot be negative.")
  
  if(!(genotype %in% x$species$genotype))
    cli::cli_abort("Genotype does not exist in the simulation.")
  
  if(!(event %in% c('death', 'growth', 'switch')))
    cli::cli_abort("Event not recognised: use any of 'death', 'growth' or 'switch'.")
  
  cc = x$simulation$get_firings() %>% 
    dplyr::filter(genotype == !!genotype, epistate == !!epistate, event == !!event) %>% 
    dplyr::pull(fired)
  
  cli::cli_alert_info("{.field {genotype}} with epistate {.field {epistate}} has count {.field {cc}} for eventt {.field {event}}, target is {.field {count}}.")
  
  x$simulation$run_up_to_event(event, paste0(genotype, epistate), count)

  return(x)
}
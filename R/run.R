#' Title
#'
#' @param x 
#' @param what 
#'
#' @return
#' @export
#'
#' @examples
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
    if('size' %in% names(what)) x = run_aux_size(x, what[['genotype']], what[['epistate']], what[['size']])
    
    if('event' %in% names(what)) x = run_aux_event(x, event = what[['event']], what[['genotype']], what[['epistate']], what[['count']])
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
    filter(genotype == !!genotype, epistate == !!epistate) %>% 
    pull(counts)
  
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
    filter(genotype == !!genotype, epistate == !!epistate, events == event) %>% 
    pull(firings)
  
  cli::cli_alert_info("{.field {genotype}} with epistate {.field {epistate}} has count {.field {cc}} for eventt {.field {event}}, target is {.field {count}}.")
  
  x$simulation$run_up_to_event(event, paste(genotype, epistate), count)
  
  return(x)
}
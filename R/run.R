run = function(
    x,
    time)
{
  stopifnot(inherits(x, 'rraces'))
  
  if(!x$has_species)
    cli::cli_abort("Species have to be added!")

    if(!x$has_initial_cell)
    cli::cli_abort("An initial cell has not been added!")
  
  if(is.null(time))
    cli::cli_abort("NULL parameters!")
  
  if(!is.numeric(time))
    cli::cli_abort("Time requires a numeric.")

  if(time < 0)
    cli::cli_abort("Time cannot be negative.")
  
  if(x$simulation$get_clock() >= time)
    cli::cli_abort("Time smaller than the current clock")
  
  cli::cli_alert_info("Current clock {.field {x$simulation$get_clock()}}, target {.field {time}}.")
  
  x$simulation$run_up_to(time)

  cli::cli_h3("New simulation state")
  
  print(x)
  
  x
  
}
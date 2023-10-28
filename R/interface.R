#' Create a simulation.
#' 
#' @description
#' This function creates a simulation that wraps calls to the RACES packace. Every
#' rRACES simulation needs to start with this function.
#' 
#' @param name A mnemonic simulation name, default "My RACES simulation".
#' @param output_dir RACES will output simulation data in this folder, by 
#' default "output_races". If the folder exists an error is raised.
#'
#' @return A simulation object
#' @export
#'
#' @examples
#' create_simulation()
create_simulation = function(
    name = "My RACES simulation",
    output_dir = 'output_races',
    tissue = "Some tissue",
    tissue_size = c(1000, 1000)
    )
{
  x = list()
  class(list) = 'rraces'
  
  if(dir.exists('output_dir'))
    cli::cli_abort("Output folder exists: {.field {output_dir}}")

  if(is.vector(tissue_size))
    cli::cli_abort("tissue_size has to be a vector with 2 entries: {.field {tissue_size}}")

  if(is.numeric(tissue_size[1]) | is.numeric(tissue_size[2]))
    cli::cli_abort("tissue_size has to be numeric: {.field {tissue_size}}")
  
  # create a new simulation and save simulation data in "output_dir"
  x$simulation <- new(Simulation, output_dir)
  x$name = name
  
  cli::cli_alert_info("Created a new simulation {.field {name}} with results in folder {.field {output_dir}}")
  
  # create a tissue
  x$simulation$set_tissue(tissue, tissue_size[1], tissue_size[2])
  cli::cli_alert_info("Added tissue {.field {tissue}} with size{.field {tissue_size[1]}x.field {tissue_size[2]}}")
  
  return(x)
}

#' Add the initial cell of the simulation.
#' 
#' @description
#' Creates a cell from which the simulation will start. The cell needs a species,
#' a epistate and a position to be placed into. Once this is inserted, no other
#' initial cells can be placed.
#'
#' @param x A simulation object.
#' @param species A valid species.
#' @param epistate An epigenetic state.
#' @param position A numeric `(x,y)` vector for the position of the cell.
#'
#' @return A simulation object.
#' @export
set_initial_cell = function(
    x,
    genotype,
    epistate,
    position
)
{
  stopifnot(inherits(x, 'rraces'))
  
  if(x$has_initial_cell)
    cli::cli_abort("An initial cell has already been added!")
  
  if(is.null(genotype) | is.null(epistate) | is.null(position))
    cli::cli_abort("NULL parameters!")
  
  if(!(genotype %in% x$species$genotype))
    cli::cli_abort("This genotype is not in the simulation")
  
  if(!(is.vector(position)))
    cli::cli_abort("Position requires a vector (x, y)")
  
  if(!is.numeric(position[1]) | !is.numeric(position[2]))
    cli::cli_abort("Position requires a numeric vector (x, y)")
  
  species_name = paste0(genotype, epistate)
  
  x$simulation$add_cell(
    species_name, 
    position[1], 
    position[2]
  )
  
  x$has_initial_cell = TRUE
  x$initial_cell = list(species_name, 
                        position[1], 
                        position[2])
  
  return(x)
}
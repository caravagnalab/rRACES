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
#'
#' @examples
set_initial_cell = function(
    x,
    species,
    epistate,
    position
)
{
  stopifnot(inherits(x, 'rraces'))
  
  if(x$has_initial_cell)
    cli::cli_abort("An initial cell has already been added!")
  
  if(is.null(species) | is.null(epistate) | is.null(position))
    cli::cli_abort("NULL parameters!")
  
  if(!(species %in% x$species$species))
    cli::cli_abort("This species is not in the simulation")
  
  if(!(is.vector(position)))
    cli::cli_abort("Position requires a vector (x, y)")
  
  if(!is.numeric(position[1]) | !is.numeric(position[2]))
    cli::cli_abort("Position requires a numeric vector (x, y)")
  
  genotype_name = paste0(species, epistate)
  
  x$simulation$add_cell(
    genotype_name = genotype_name, 
    x = position[1], 
    y = position[2]
  )
  
  x$has_initial_cell = TRUE
  
  return(x)
}
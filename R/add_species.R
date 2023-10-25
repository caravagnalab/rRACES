#' Add a species
#'
#' @param x A simulation object
#' @param transition_rate Transition rates from (+) to (-), and viceversa
#' @param growth_rate Growth rates from (+) to (-), and viceversa
#' @param death_rate
#'
#' @return
#' @export
#'
#' @examples
add_species = function(x,
                       transition_rate = c(`+-` = 0.01, `-+` = 0.01),
                       growth_rate = c(`+` = 0.2, `-` = 0.15),
                       death_rate = c(`+` = 0.1, `-` = 0.1)
                       )
{
  stopifnot(class(x) == 'w_rraces')

  # apply to X interface functions

}

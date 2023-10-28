#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @exportS3Method print.Simulation
#'
#' @examples
print.rraces <- function(x, ...) {
  
  cli::cli_h1(paste("rRACES:", x$name))
  
  sp_name = x$simulation$get_species_names()
  
  cli::cli_alert("Species {.field {sp_name}}")
  
  invisible(x)
}







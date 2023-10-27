#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @exportS3Method print.Simulation
#'
#' @examples
print.Simulation <- function(x, ...) {
  cat("Object with value:", x$get_genotype_names(), "\n")
  invisible(x)
}







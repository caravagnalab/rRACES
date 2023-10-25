#' Print to screen a simulation object.
#'
#' @param x simulation object.
#'
#' @exportS3Method
#'
#' @examples
#'
#' # will automatically trigger print
#' init()
print.w_rraces = function(x, ...){
  cli::cli_h1(x$name)
}

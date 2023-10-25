#' Plot a simulation object.
#'
#' @param x simulation object.
#'
#' @exportS3Method
#'
#' @examples
#'
#' # will automatically trigger print
#' plot(init())
plot.w_rraces = function(x, ...){
  ggplot2::ggplot()
}

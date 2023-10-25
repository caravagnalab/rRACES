#' Construct a simulation object, with a name.
#'
#' @param name simulation name.
#'
#' @return an S3 object of class w_rraces
#' @export
#'
#' @examples
#'
#' init()
#' init("example simul.")
init = function(name = "My RACES simulation"){

  x = list()
  class(x) = 'w_rraces'

  x$name = name

  return(x)
}

#' Create a tissue.
#'
#' @param name The name of the tissue
#' @param width Width in pixels
#' @param height Height in pixels
#' @param x A simulation object
#'
#' @return
#' @export
#'
#' @examples
#' init() %>% create_tissue()
create_tissue = function(x, name = "My Tissue", width = 1000, height = 1000)
{
  stopifnot(class(x) == 'w_rraces')

  # apply to X interface functions

}

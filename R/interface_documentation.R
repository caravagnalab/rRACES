#' #' Title
#' #'
#' #' @param simulation 
#' #' @param tissue 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' set_tissue <- function(simulation, tissue) {
#'   .Call("Simulation_set_tissue", simulation, tissue)
#' }
#' 
#' #' Title
#' #'
#' #' @param simulation 
#' #' @param species 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' add_species <- function(simulation, species) {
#'   .Call("Simulation_add_species", simulation, species)
#' }
#' 
#' #' Title
#' #'
#' #' @param simulation 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' get_genotype_names <- function(simulation) {
#'   .Call("Simulation_get_genotype_names", simulation)
#' }
#' 
#' #' Title
#' #'
#' #' @param simulation 
#' #' @param cell 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' add_cell <- function(simulation, cell) {
#'   .Call("Simulation_add_cell", simulation, cell)
#' }
#' 
#' #' Title
#' #'
#' #' @param simulation 
#' #' @param time 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' run_up_to <- function(simulation, time) {
#'   .Call("Simulation_run_up_to", simulation, time)
#' }
#' Add genotype evolution.
#' 
#' @description
#' This function adds a genotype evolution event to a simulation object. Genotype 
#' evolution is the transition from one genotype to another at a specified time. 
#' The function accepts the simulation object and the start/end genotypes, and the
#' `time` at which the event happens. When the clock reaches `time`, a cell is
#' randomly sampled from those with the starting genotype (`genotype_from`) and
#' its genotype is changed to `genotype_to`. Note that the epistate of the 
#' cell is therefore sampled at randome based on the abundance of the cells.
# 
#' @param x An 'rraces' simulation object.
#' @param genotype_from The starting genotype before the evolution.
#' @param genotype_to The target genotype after the evolution.
#' @param time The time at which the genotype evolution occurs.
#'
#' @return An updated simulation object with the added genotype evolution event.
#'
#' @export
#'
#' @examples
#'
#' # Example of using the function
#' x <- create_simulation()
#' # Add two species, "A" and "B", all default parameters
#' x = add_species(x, "A")
#' x = add_species(x, "B")
#' #Program a timed transition 
#' x = add_genotype_evolution(x, "A", "B", 10)
#' print(x)
add_genotype_evolution = function(x,
                                  genotype_from,
                                  genotype_to,
                                  time)
{
  stopifnot(inherits(x, 'rraces'))
  
  
  # Check inputs
  if (!x$has_species)
    cli::cli_abort("Need to add species first")
  
  if (is.null(genotype_from) | is.null(genotype_to))
    cli::cli_abort("Genotypes cannot be NULL")
  
  if (!(genotype_from %in% x$species$genotype))
    cli::cli_abort("Genotype {.value {genotype_from}} does not exists.")
  
  if (!(genotype_to %in% x$species$genotype))
    cli::cli_abort("Genotype {.value {genotype_to}} does not exists.")
  
  t_ancestor = TimeOf(x, genotype_from)
  
  if (!is.na(t_ancestor) & t_ancestor >= time)
    cli::cli_abort(
      "{.value {genotype_from}} is born at time {.value {t_ancestor}}, {.value {genotype_to}} cannot be born at time {.field {time}}"
    )

  # We can add it now
  x$simulation$schedule_genotype_mutation(genotype_from, genotype_to, time)
  
  x$species = x$species %>% 
    dplyr::mutate(
      ancestor = ifelse(
        genotype == genotype_to,
        genotype_from,
        ancestor
      ),
      time = ifelse(
        genotype == genotype_to,
        !!time,
        time
      )
    )
  
  df = x$species %>% 
    dplyr::select(ancestor, genotype, time) %>% 
    dplyr::distinct() %>% 
    dplyr::rename(from = ancestor, to = genotype)
  
  which_roots = df %>% dplyr::filter(is.na(from)) %>% pull(to)
  
  df = rbind(df,
                 data.frame(
                   from = 'GL',
                   to = which_roots,
                   time = 0,
                   stringsAsFactors = FALSE
                 )) %>% 
    dplyr::filter(!is.na(from))
  
  x$clone_tree = tidygraph::as_tbl_graph(df)
  
  x$has_species = TRUE

  cli::cli_alert_info(
    "{crayon::underline('New transition')} {.field {genotype_from}} -> {.field {genotype_to}} at time {.field {time}}"
  )

  print(x)
}
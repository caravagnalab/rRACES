#' Title
#'
#' @param x 
#' @param genotype_from 
#' @param genotype_to 
#' @param time 
#'
#' @return
#' @export
#'
#' @examples
add_genotype_evolution = function(x,
                                  genotype_from,
                                  genotype_to,
                                  time)
{
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
  x$simulation$add_timed_mutation(genotype_from, genotype_to, time)
  
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
    dplyr::distinct %>% 
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
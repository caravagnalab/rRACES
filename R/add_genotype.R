#' Add a genotype to a simulation.
#'
#' @description
#' This function adds one new genotype to the simulation.
#' 
#' @param x A simulation.
#' @param name The genotype name.
#' @param epigenetic_rates The epigenetic rates of the new species.
#' @param growth_rates The duplication rates of the new species.
#' @param death_rates The death rates of the new species.
#'
#' @return A simulation object.
#' @export
#'
#' @examples
#' x <- create_simulation(output_dir="add_genotype_test")
add_genotype = function(x,
                        name = "A",
                        epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
                        growth_rates = c("+" = 0.2, "-" = 0.08),
                        death_rates = c("+" = 0.1, "-" = 0.01))
{
  stopifnot(inherits(x, 'rraces'))
  
  # Check inputs
  if (name %in% x$species$genotype)
    cli::cli_abort("Species {.value {name}} already exists.")
  
  if (is.null(epigenetic_rates))
    cli::cli_abort("Epigenetic rates are NULL")
  else{
    if (!all(names(epigenetic_rates) %in% c("+-", "-+")))
      cli::cli_abort("Epigenetic rates are missing fields")
    if (any(epigenetic_rates < 0))
      cli::cli_abort("Epigenetic rates are negative")
  }
  
  if (is.null(growth_rates))
    cli::cli_abort("Growth rates are NULL")
  else{
    if (!all(names(growth_rates) %in% c("+", "-")))
      cli::cli_abort("Growth rates are missing fields")
    if (any(growth_rates < 0))
      cli::cli_abort("Growth rates are negative")
  }
  
  if (is.null(death_rates))
    cli::cli_abort("Death rates are NULL")
  else{
    if (!all(names(death_rates) %in% c("+", "-")))
      cli::cli_abort("Death rates are missing fields")
    if (any(death_rates < 0))
      cli::cli_abort("Death rates are negative")
  }
  
  # We can add it now
  x$simulation$add_genotype(
    name = name,
    epigenetic_rates = epigenetic_rates,
    growth_rates = growth_rates,
    death_rates = death_rates
  )
  
  x$species = x$species %>%
    dplyr::bind_rows(
      dplyr::tibble(
        species = paste0(name, '+'),
        genotype = name,
        ancestor = NA,
        time = NA,
        epistate = '+',
        rgrowth = growth_rates['+'],
        rdeath = death_rates['+'],
        repigenetic = epigenetic_rates['+-']
      )
    ) %>%
    dplyr::bind_rows(
      dplyr::tibble(
        species = paste0(name, '-'),
        genotype = name,
        ancestor = NA,
        time = NA,
        epistate = '-',
        rgrowth = growth_rates['-'],
        rdeath = death_rates['-'],
        repigenetic = epigenetic_rates['-+']
      )
    )
  
  x$has_species = TRUE
  
  cli::cli_alert_info(
    "{crayon::underline('New species')} {.field {name}}, added."
  )
  
  print(x)
  
  return(x)
}
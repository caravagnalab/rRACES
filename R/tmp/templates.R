
#' Simulation: template monoclonal tumour without epigenetic events
#'
#' @return A JSON valid input for RACES
#' @export
#'
#' @examples
#' template_monoclonal_noepi()
template_monoclonal_noepi = function()
{
  # Monoclonal A- population
  genotypes = list(
    list(
      name = "A",
      "epigenetic rates" = list(
        list(
          methylation = 0,
          demethylation = 0
        )
      ),
      "epigenetic types" = list(
        list(
          status = "-",
          "death rate" = 0.1,
          "duplication rate" = 0.2
        )
      )
    )
  )

  initial_cells <- list(
    list(
      genotype = list(
        name = "A",
        "epigenetic status" = "-"
      ),
      position = c(500, 500)
    )
  )

  timed_events <- list(
    list(
      time = 100,
      type = "sampling",
      "remove sample" = TRUE,
      "sample regions" = list(
        list("lower corner" = c(450, 450), "upper corner" = c(455, 455)),
        list("lower corner" = c(540, 540), "upper corner" = c(545, 545))
      )
    )
  )

  # Call the function to generate JSON
  json_data <- generate_custom_json(
    comment = "Template monoclonal tumour without epigenetic events",
    genotypes = genotypes,
    tissue_name = "Custom Tissue",
    tissue_size = c(1200, 1200),
    initial_cells = initial_cells,
    timed_events = timed_events,
    death_activation_level = 150
  )

  return(json_data)
}

#' Simulation: template monoclonal tumour with epigenetic events
#'
#' @return A JSON valid input for RACES
#' @export
#'
#' @examples
#' template_monoclonal_epi()
template_monoclonal_epi = function()
{
  # Monoclonal A-/+ population
  genotypes = list(
    list(
      name = "A",
      "epigenetic rates" = list(
        list(
          methylation = 0.01,
          demethylation = 0.01
        )
      ),
      "epigenetic types" = list(
        list(
          status = "-",
          "death rate" = 0.1,
          "duplication rate" = 0.2
        ),
        list(
          status = "+",
          "death rate" = 0.05,
          "duplication rate" = 0.2
        )
      )
    )
  )

  initial_cells <- list(
    list(
      genotype = list(
        name = "A",
        "epigenetic status" = "-"
      ),
      position = c(500, 500)
    )
  )

  timed_events <- list(
    list(
      time = 100,
      type = "sampling",
      "remove sample" = TRUE,
      "sample regions" = list(
        list("lower corner" = c(450, 450), "upper corner" = c(455, 455)),
        list("lower corner" = c(540, 540), "upper corner" = c(545, 545))
      )
    )
  )

  # Call the function to generate JSON
  json_data <- generate_custom_json(
    comment = "Template monoclonal tumour with epigenetic events",
    genotypes = genotypes,
    tissue_name = "Custom Tissue",
    tissue_size = c(1200, 1200),
    initial_cells = initial_cells,
    timed_events = timed_events,
    death_activation_level = 150
  )

  return(json_data)
}

#' Simulation: template tumour 2 clones and no epigenetic events
#'
#' @return A JSON valid input for RACES
#' @export
#'
#' @examples
#' template_2clones_noepi()
template_2clones_noepi = function()
{
  # Example usage:
  genotypes <- list(
    list(
      name = "A",
      "epigenetic rates" = list(
        list(
          methylation = 0,
          demethylation = 0
        )
      ),
      "epigenetic types" = list(
        list(
          status = "-",
          "death rate" = 0.1,
          "duplication rate" = 0.2
        )
      )
    ),
    list(
      name = "B",
      "epigenetic rates" = list(
        list(
          methylation = 0,
          demethylation = 0
        )
      ),
      "epigenetic types" = list(
        list(
          status = "-",
          "death rate" = 0.1,
          "duplication rate" = 0.3
        )
      )
    )
  )

  initial_cells <- list(
    list(
      genotype = list(
        name = "A",
        "epigenetic status" = "-"
      ),
      position = c(500, 500)
    )
  )

  timed_events <- list(
    list(
      time = 100,
      type = "sampling",
      "remove sample" = TRUE,
      "sample regions" = list(
        list("lower corner" = c(450, 450), "upper corner" = c(455, 455)),
        list("lower corner" = c(540, 540), "upper corner" = c(545, 545))
      )
    ),
    list(
      time = 60,
      type = "driver mutation",
      "original genotype" = "A",
      "mutated genotype" = "B"
    )
  )

  json_data <- generate_custom_json(
    comment = "Template 2-clones tumour without epigenetic events",
    genotypes = genotypes,
    tissue_name = "Custom Tissue",
    tissue_size = c(1200, 1200),
    initial_cells = initial_cells,
    timed_events = timed_events,
    death_activation_level = 150
  )

  return(json_data)
}

#' Simulation: template tumour 2 clones with epigenetic events
#'
#' @return A JSON valid input for RACES
#' @export
#'
#' @examples
#' template_2clones_epi()
template_2clones_epi = function()
{
  # Example usage:
  genotypes <- list(
    list(
      name = "A",
      "epigenetic rates" = list(
        list(
          methylation = 0.01,
          demethylation = 0.01
        )
      ),
      "epigenetic types" = list(
        list(
          status = "-",
          "death rate" = 0.1,
          "duplication rate" = 0.2
        ),
        list(
          status = "+",
          "death rate" = 0.05,
          "duplication rate" = 0.2
        )
      )
    ),
    list(
      name = "B",
      "epigenetic rates" = list(
        list(
          methylation = 0.01,
          demethylation = 0.01
        )
      ),
      "epigenetic types" = list(
        list(
          status = "-",
          "death rate" = 0.1,
          "duplication rate" = 0.3
        ),
        list(
          status = "+",
          "death rate" = 0.05,
          "duplication rate" = 0.3
        )
      )
    )
  )

  initial_cells <- list(
    list(
      genotype = list(
        name = "A",
        "epigenetic status" = "-"
      ),
      position = c(500, 500)
    )
  )

  timed_events <- list(
    list(
      time = 100,
      type = "sampling",
      "remove sample" = TRUE,
      "sample regions" = list(
        list("lower corner" = c(450, 450), "upper corner" = c(455, 455)),
        list("lower corner" = c(540, 540), "upper corner" = c(545, 545))
      )
    ),
    list(
      time = 60,
      type = "driver mutation",
      "original genotype" = "A",
      "mutated genotype" = "B"
    )
  )

  json_data <- generate_custom_json(
    comment = "Template 2-clones tumour with epigenetic events",
    genotypes = genotypes,
    tissue_name = "Custom Tissue",
    tissue_size = c(1200, 1200),
    initial_cells = initial_cells,
    timed_events = timed_events,
    death_activation_level = 150
  )

  return(json_data)
}

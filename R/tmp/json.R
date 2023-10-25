generate_custom_json <- function(
    comment = "An input example for \"drivers_sim\".",
    genotypes,
    tissue_name = "Liver",
    tissue_size = c(1000, 1000),
    initial_cells,
    timed_events,
    death_activation_level = 100
) {
  data <- list(
    comment = comment,
    genotypes = genotypes,
    tissue = list(
      name = tissue_name,
      size = tissue_size
    ),
    "initial cells" = initial_cells,
    "timed events" = timed_events,
    "death activation level" = death_activation_level
  )

  # Convert the data to JSON
  json_data <- jsonlite::toJSON(data, pretty = TRUE)

  return(json_data)
}

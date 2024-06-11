
#' Get stick table
#'
#' @description
#' Get subset of forest nodes corresponding to phylogenetic branches containing
#'  relevant biological events (drivers)
#'
#' @param forest The original forest object
#' has been derived
#'
#' @return A `tidyverse` tibble object
#' @export
#'
#' @examples
#' sim <- new(Simulation)
#'
#' sim$add_mutant(name = "A",
#'                growth_rates = 1,
#'                death_rates = 0)
#' sim$place_cell("A", 500, 500)
#'
#' sim$run_up_to_size("A",1e4)
#'
#' sim$add_mutant(name = "B",
#'                growth_rates = 3.5,
#'                death_rates = 0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#'
#' sim$run_up_to_size("B",1e4)
#'
#' bbox <- sim$search_sample(c("A" = 100,"B" = 100), 50, 50)
#' sim$sample_cells("Sampling", bbox$lower_corner, bbox$upper_corner)
#'
#' forest = sim$get_samples_forest()
#' get_events_table(forest)

get_events_table = function(forest) {

  sticks <- forest$get_sticks()

  if (length(sticks) > 0) {
    events_table <- lapply(1:length(sticks), function(i) {

      gen <- forest$get_nodes() %>% dplyr::as_tibble() %>%
        filter(cell_id == sticks[[i]][length(sticks[[i]])]) %>% pull(mutant)

      epi <- forest$get_nodes() %>% dplyr::as_tibble() %>%
        filter(cell_id == sticks[[i]][length(sticks[[i]])]) %>% pull(epistate)

      if (epi != "") {
        change <- paste0(ifelse(epi == "+", "-", "+"), " -> ", epi)
        dplyr::tibble(cell_id = sort(sticks[[i]][2:length(sticks[[i]])]),
                      label = paste0(gen, " ", change,
                                     "(",
                                     sort(sticks[[i]])[length(sticks[[i]])],
                                     ")"))
      } else{
        dplyr::tibble(cell_id = sort(sticks[[i]][2:length(sticks[[i]])]),
                      label = gen)
      }

    }) %>% dplyr::bind_rows()

    dplyr::as_tibble(forest$get_nodes()) %>%
      dplyr::full_join(events_table, by = "cell_id") %>%
      mutate(label = ifelse(is.na(label),"Subclonal",label)) %>%
      mutate(label = ifelse(cell_id == 0,"Truncal",label)) %>%
      unique()

  }else{
    NULL
  }

}


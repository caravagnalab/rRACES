## This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
## Copyright (C) 2023 - Alberto Casagrande <alberto.casagrande@uniud.it>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' @importFrom rlang .data
setMethod("show", "Rcpp_Simulation", function(object) {
  # If it can be simulated
  sim_status <- (length(object$get_species())
                 & nrow(object$get_cells()))
  
  # If it has samples assigned
  sam_status <- sim_status & nrow(object$get_samples_info())
  
  mut_status <- FALSE

  sim_status <- ifelse(
    sim_status,
    crayon::bgGreen(crayon::white(" D ")),
    crayon::bgRed(crayon::white(" D "))
  )

  sam_status <- ifelse(
    sam_status,
    crayon::bgGreen(crayon::white(" S ")),
    crayon::bgRed(crayon::white(" S "))
  )

  mut_status <- ifelse(
    mut_status,
    crayon::bgGreen(crayon::white(" M ")),
    crayon::bgRed(crayon::white(" M "))
  )

  cli::cli_rule(
    left = paste(
      crayon::bgYellow(" rRACES "),
      sim_status, sam_status, mut_status,
      object$get_name()
    ),
    right = paste0(
      crayon::red("\u25A3 "),
      object$get_tissue_name(),
      " [",
      paste(object$get_tissue_size(), collapse = "x"),
      "]",
      crayon::red("  \u23f1 "),
      signif(object$get_clock(), digits = 3)
    )
  )

  species <- object$get_species()
  if (nrow(species) == 0) {
    cat('\n')
    cli::cli_alert_danger("The simulation has no species yet!")
  } else {
    counts_tab <- object$get_counts() %>% 
      dplyr::mutate(
        `%` = 100 * counts/sum(counts)
      )

    my_tab <- species %>%
      dplyr::left_join(counts_tab, by = c("genotype", "epistate")) %>%
      dplyr::mutate(species = paste0(.data$genotype, .data$epistate)) %>%
      dplyr::select(.data$species, .data$growth_rate,
                    .data$death_rate, .data$switch_rate, .data$counts, .data$`%`)
    colnames(my_tab)[2:4] <- c(" \u03BB ", " \u03B4 ", " \u03B5 ")
    
    nspecies = nrow(my_tab)
    with_epigenetics = max(nchar(species$epistate)) != 0
    
    cli::cli_h3(text = "Species: {.field {nspecies}}, {.field {ifelse(with_epigenetics, crayon::green('with'), crayon::red('without'))}} epigenetics")

    cat("   ",
        knitr::kable(my_tab, format = "rst", align = "rcrccc"),
        sep = "\n   ")
    
    # cli::cli_rule(center = "Firings")
    
    f_table = object$get_firings() %>%
      dplyr::mutate(species = paste0(.data$genotype, .data$epistate))
    
    cli::cli_h3(text = paste("Firings:", sum(f_table$fired), 'total'))
    
    nch = f_table$fired %>% nchar %>% max
    
    for(s in my_tab$species)
    {
      s_ftab = f_table %>% dplyr::filter(species == s) %>% dplyr::arrange(event)
      
      if(with_epigenetics)
        sprintf(
          paste0('\n\tSpecies [%s]: %', nch, 's (death), %', nch, 's (growth) and %', nch, 's (switches)'),
          s,
          s_ftab$fired[1],
          s_ftab$fired[2],
          s_ftab$fired[3]
        ) %>% cat()
      else
          sprintf(
            paste0('\n\tSpecies [%s]: %', nch, 's (death) and %', nch, 's (growth)'),
            s,
            s_ftab$fired[1],
            s_ftab$fired[2]
          ) %>% cat()
      
    }
    
    # cat("   ",
    #     knitr::kable(object$get_firings(), format = "rst", align = "rcrccc"),
    #     sep = "\n   ")
    
  }
  
  
  samples <- object$get_samples_info()
  if (nrow(samples) == 0) {
    cat('\n')
    cli::cli_alert_danger("The simulation has no samples yet!")
  } else {
    t_samples = samples$time %>% unique()
    
    is_longitudinal = length(t_samples) > 1
    is_multiregion = any(table(samples$time) > 1)
    
    is_longitudinal = ifelse(is_longitudinal, crayon::bgGreen(crayon::white(" Yes ")), crayon::bgRed(crayon::black(" No ")))
    is_multiregion = ifelse(is_multiregion, crayon::bgGreen(crayon::white(" Yes ")), crayon::bgRed(crayon::black(" No ")))

    cat("\n")
    
    ncells = samples$tumoural.cells %>% sum
    
    cli::cli_h3(text = "Samples ({.emph multi-region} {is_multiregion}, {.emph longitudinal} {is_longitudinal}): {.field {ncells}} cells from {.field {nrow(samples)}} samples at {.field {length(t_samples)}} timepoints")
    # 
    # 
    # cli::cli_alert_info(
    #   "There are {.field {nrow(samples)}} samples from {.field {length(t_samples)}} timepoints."
    # )
    
    cat("   ",
        knitr::kable(samples, format = "rst", align = "rcrccc"),
        sep = "\n   ")
 
  }
})
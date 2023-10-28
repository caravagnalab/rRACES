#' S3 print object for rRACES.
#'
#' @param x An object returned by `create_simulation`.
#' @param ... Ellipsis.
#'
#' @return
#' @exportS3Method print.rraces
#'
#' @examples
print.rraces <- function(x, ...) {
  
  cli::cli_rule(
    left = paste(crayon::bold("rRACES:"), x$name),
    righ = paste("⏱ ", x$simulation$get_clock())
  )
  cat('\n')
  
  # Tissue
  # 
  # x$simulation$get_clock()
  cli::cli_alert_warning("TODO -- print the tissue with get_info")
  
  # Species
  sp_name = x$species$species %>% unique()
  if(!x$has_species)
    cli::cli_alert_danger("The simulation has no species yet!")
  else
  {
    # Pretty print
    max_pad = 1 + (nchar(x$species$species) %>% max)
    max_pad = ifelse(max_pad <=2, 3, max_pad)
    max_pad = paste0('%', max_pad, 's')
    
    for(s in sp_name)
    {
      cnt = x$simulation$get_counts() %>% filter(species == s)
      
      sname = paste('', s, '')
      sname = sprintf(max_pad, s)
      sname = crayon::bgYellow(crayon::black( sname ))
      
      ancestor = AncestorOf(x, s)
      ancestor = ifelse(is.na(ancestor), 'wt', ancestor)
      ancestor = sprintf(max_pad, ancestor)
      ancestor = crayon::bgBlue(crayon::white( ancestor ))
      ancestor = paste0('\u21AA ', ancestor)
      
      grate = paste0(
        '  \u03BB =', sprintf("%6s", GRateOf(x, s, '+')), ' (+) ', sprintf("%6s", GRateOf(x, s, '-')), ' (-)'
      )

      drate = paste0(
        '  \u03B4 =', sprintf("%6s", DRateOf(x, s, '+')), ' (+) ', sprintf("%6s", DRateOf(x, s, '-')), ' (-)'
      )

      erate = paste0(
        '  \u03B5 =', sprintf("%6s", ERateOf(x, s, '+')), ' (+) ', sprintf("%6s", ERateOf(x, s, '-')), ' (-)'
      )
      
      cat('\n\t\u2022', sname, ancestor, "|", grate, "|", drate, "|", erate)   
    }
  }
  
  can_simulate = FALSE
  can_generate_data = FALSE
  
  invisible(x)
}







#' S3 print for rRACES object.
#' @description
#' S3 print method.
#'
#' @param x An object returned by `create_simulation`.
#' @param ... Ellipsis.
#'
#' @export
#' 
#' @examples
#' 
#' # Example of using the function
#' x = create_simulation()
#' x = add_species(x, "A")
#' x = set_initial_cell(x, "A", "+", c(50, 50))
#' x = run(x, list(time = 60))
#' print(x)
print.rraces <- function(x, ...) {
  
  SIM_status = (x$has_species & x$has_initial_cell)
  SAM_status = FALSE
  MUT_status = FALSE
  
  SIM_status = ifelse(
    SIM_status,
    crayon::bgGreen(crayon::white(' D ')),
    crayon::bgRed(crayon::white(' D '))
  )  
  
  SAM_status = ifelse(
    SAM_status,
    crayon::bgGreen(crayon::white(' S ')),
    crayon::bgRed(crayon::white(' S '))
  )  
  
  MUT_status = ifelse(
    MUT_status,
    crayon::bgGreen(crayon::white(' M ')),
    crayon::bgRed(crayon::white(' M '))
  )  
  
  cli::cli_rule(
    left = paste(
      crayon::bgYellow(" rRACES "), 
      SIM_status, SAM_status, MUT_status,
      x$name
      ),
    righ = paste0(
      crayon::red("\u25A3 "),
      x$simulation$get_tissue_name(),
      " [",
      paste(x$simulation$get_tissue_size(), collapse = 'x'),
      "]",
      crayon::red("\t\u23F1 "),
      x$simulation$get_clock()
    )
  )
  # cat('\n')
  
  # Species
  sp_name = x$species$species %>% unique()
  
  if (!x$has_species)
  {
    cat('\n')
    cli::cli_alert_danger("The simulation has no species yet!")
  }
  else
  {
    # Pretty print
    # max_pad = 1 + (nchar(x$species$genotype) %>% max)
    # max_pad = ifelse(max_pad <=2, 3, max_pad)
    # max_pad = paste0('%', max_pad, 's')
    #
    # for(s in sp_name)
    # {
    #   cnt = x$simulation$get_counts() %>% filter(genotype == s)
    #
    #   sname = paste('', s, '')
    #   sname = sprintf(max_pad, s)
    #   sname = crayon::bgYellow(crayon::black( sname ))
    #
    #   ancestor = AncestorOf(x, s)
    #   ancestor = ifelse(is.na(ancestor), 'wt', ancestor)
    #   ancestor = sprintf(max_pad, ancestor)
    #   ancestor = crayon::bgBlue(crayon::white( ancestor ))
    #   ancestor = paste0('\u21AA ', ancestor)
    #
    #   grate = paste0(
    #     '  \u03BB =', sprintf("%6s", GRateOf(x, s, '+')), ' (+) ', sprintf("%6s", GRateOf(x, s, '-')), ' (-)'
    #   )
    #
    #   drate = paste0(
    #     '  \u03B4 =', sprintf("%6s", DRateOf(x, s, '+')), ' (+) ', sprintf("%6s", DRateOf(x, s, '-')), ' (-)'
    #   )
    #
    #   erate = paste0(
    #     '  \u03B5 =', sprintf("%6s", ERateOf(x, s, '+')), ' (+) ', sprintf("%6s", ERateOf(x, s, '-')), ' (-)'
    #   )
    #
    #   cat('\n\t\u2022', sname, ancestor, "|", grate, "|", drate, "|", erate)
    # }
    #
    # library(knitr)
    #
    
    # Call the function with your tibble 'x$species'
    counts_tab = x$simulation$get_counts()
    
    my_tab = x$species %>% 
      dplyr::left_join(counts_tab, by = c('genotype', 'epistate')) %>% 
      dplyr::select(species, ancestor, time, rgrowth, rdeath, repigenetic, counts)
    colnames(my_tab)[4:6] = c(" \u03BB ", " \u03B4 ", " \u03B5 ")
    
    cat("   ",
        knitr::kable(my_tab, format = "rst", align = "rcrccc"),
        sep = "\n   ")
  }
  
  if(x$has_initial_cell)
  {
    cat(
      "\n  Initial cell:",
      crayon::blue(x$initial_cell[1]),
      '→',
      paste0(
        '(',
        crayon::blue(x$initial_cell[2]),
        'x',
        crayon::blue(x$initial_cell[3]),
        ')'
      ),
      '\n'
    )
  }
  
  
  # cat("\n")
  # if (x$has_species & x$has_initial_cell)
  #   cat(crayon::bgGreen(crayon::white(' RACES can simulate ')))
  # else
  #   cat(crayon::bgRed(crayon::white(' RACES cannot simulate ')))
  
  invisible(x)
}


# printPretty = function(node, indent, last)
# {
#   cat(indent)
#   if (last) {
#     cat("\\-")
#     indent = paste(indent, " ", sep = '')
#   }
#   else {
#     cat("|-")
#     indent = paste(indent, "| ", sep = '')
#   }
#   cat(node)
#   
#   A = x$clone_tree %>%
#     tidygraph::activate(nodes) %>%
#     dplyr::filter(name == !!node) %>%
#     pull(time)
# 
#   
#   cat('\n')
#   
#   cl = children(M, node)
#   
#   for (c in cl)
#     printPretty(c, indent, c == cl[length(cl)])
# }


# Plot call
# layout <- ggraph::create_layout(x$clone_tree, layout = 'tree')
# 
# ggraph(layout) +
#   geom_edge_link(
#     arrow = arrow(length = unit(2, 'mm')),
#     end_cap = circle(5, 'mm'),
#     start_cap  = circle(5, 'mm')
#   ) +
#   geom_node_point(aes(colour = name,
#                       size = time),
#                   na.rm = TRUE) +
#   geom_node_text(aes(label = name),
#                  colour = 'black',
#                  vjust = 0.4) +
#   coord_cartesian(clip = 'off') +
#   # theme_graph(base_size = 8, base_family = '') +
#   theme_void(base_size = 10) +
#   theme(legend.position = 'bottom') +
#   scale_color_manual(values = get_species_colors(x$species$genotype %>% unique)) +
#   guides(color = FALSE,
#          size = guide_legend("Clone size", nrow = 1)) 


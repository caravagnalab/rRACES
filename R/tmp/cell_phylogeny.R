# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#
# ## BiocManager::install("BiocUpgrade") ## you may need this
# BiocManager::install("treeio")

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
parse_cell_phylogeny = function(
    x = '~/Downloads/RACES/phylo.xml'
      )
{
  tree = treeio::read.phyloxml('~/Downloads/RACES/phylo.xml')

  library(ggtree)
  library(ggplot2)
  library(dplyr)

  tree@data = tree@data %>%
    mutate(
      col = paste(red, green, blue)
    )

  ggtree(tree) +
    geom_tiplab(size = 1.5)  + # Add labels to the tips (leaves)
    geom_nodepoint(aes(color = col), size = 1) +
    scale_color_brewer(palette = 'Set2') +
    theme(
      legend.position = 'bottom'
    )

}



# tree = x
# ggtree(tree, layout="slanted")
# ggtree(tree, layout="circular")
# ggtree(tree, layout="fan", open.angle=120)
# ggtree(tree, layout="equal_angle")
# ggtree(tree, layout="daylight")
# ggtree(tree, branch.length='none')
# ggtree(tree, branch.length='none', layout='circular')
# ggtree(tree, layout="daylight", branch.length = 'none')
# ggtree(tree, mrsd="2013-01-01") + theme_tree2()
#
# ggtree(tree) + geom_tiplab() +
#   geom_label(aes(x=branch, label=red), fill='lightgreen') +
#   geom_label(aes(label=green), fill='steelblue') +
#   geom_text(aes(label=blue), hjust=-.5)


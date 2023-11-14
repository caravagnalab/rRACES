library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

# Set the number of vertices in the tree
num_vertices <- 1500

# Generate a random tree
random_tree <- sample_tree(num_vertices,directed = TRUE)

# Get the list of edges
edge_list <- get.data.frame(random_tree, what = "edges")

isLeaf = function(x, edge_list){
  (edge_list %>% dplyr:: filter(from == x) %>% nrow()) == 0
}

incomeEdges = function(x, edge_list){
  edge_list %>% dplyr:: filter(to == x)
}

nodes = function(edge_list){
  edge_list %>% unlist %>% unique 
}

RootNode = function(edge_list)
{
  istar = sapply(nodes(edge_list), function(w){
    incomeEdges(w, edge_list) %>% nrow()
  })
  names(istar) = nodes
  names(which(istar == 0))
}


nodes_list = data.frame(
  nodeID = nodes(edge_list),
  stringsAsFactors = FALSE
) %>% 
  rowwise() %>% 
  dplyr::mutate(
    sample =  ifelse(isLeaf(nodeID, edge_list), "BBOX", NA),
    genotype = sample(c("A", "B"), 1),
    epistate = sample(c("+", "-"), 1)
  )

# Create a tidygraph object
graph <- as_tbl_graph(edge_list, directed = TRUE)

# Activate nodes and add node attributes
graph <- graph %>%
  activate(nodes) %>%
  mutate(
    sample = nodes_list$sample,
    genotype = nodes_list$genotype,
    epistate = nodes_list$epistate,
    species = paste0(genotype, epistate)
  )

# Define the layout with the tree rooted and the root on top
layout <- create_layout(graph, layout = "tree", root = 1)

# species_colors = get_species_colors(x)
species_colors = RColorBrewer::brewer.pal(4, "Set1")

ncells = graph %>%
  activate(nodes) %>% 
  pull(name) %>% 
  length

nsamples = graph %>%
  activate(nodes) %>% 
  filter(!is.na(sample)) %>% 
  pull(name) %>% 
  length

labels_every = (layout$y %>% max)/10

# Plot the graph with color-coded nodes
ggraph(layout, "tree") +
  geom_edge_link() +
  geom_node_point(
    aes(
      color = species, 
      shape = ifelse(is.na(sample), "N/A", sample)
      ), 
      size = 2) +
  scale_color_manual(values = species_colors) +
  # scale_shape_manual(values = c("+" = 16, "-" = 17)) +
  # theme_void() +
  theme_minimal() + 
  theme(legend.position = "bottom") +
  labs(
    title = "Cell divison tree",
    subtitle = paste0(ncells, " cells, ", nsamples, " samples"),
    shape = "Sample",
    x = NULL,
    y = "Cell division"
  ) +
  scale_y_continuous(labels = 
                       seq(0, layout$y %>% max, labels_every) %>% round %>% rev,
                     breaks = 
                       seq(0, layout$y %>% max, labels_every) %>% round
                     ) +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )





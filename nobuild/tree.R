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




library(rRACES)
library(dplyr)
sim <- new(Simulation, "Monoclonal")

sim$add_genotype(name = "A",
                 epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
                 growth_rates = c("+" = 0.1, "-" = 0.08),
                 death_rates = c("+" = 0.1, "-" = 0.01))

sim$place_cell("A+", 500, 500)
sim$run_up_to_size("A+", 1000)

plot_tissue(sim, num_of_bins = 500)

# a 5x5 box gives 25 cells
bbox_width = 10

sim$sample_cells("A1", bottom_left = c(480, 480), top_right = c(480 + bbox_width, 480 + bbox_width))
sim$sample_cells("B1", bottom_left = c(500, 500), top_right = c(500 + bbox_width, 500 + bbox_width))

plot_tissue(sim, num_of_bins = 500)

sim$run_up_to_time(sim$get_clock() + 15)

plot_tissue(sim, num_of_bins = 500)

cell <- sim$choose_cell_in("A")

sim$add_genotype(name = "B",
                 epigenetic_rates = c("+-" = 0.05, "-+" = 0.1),
                 growth_rates = c("+" = 0.5, "-" = 0.3),
                 death_rates = c("+" = 0.05, "-" = 0.05))

sim$mutate_progeny(cell, "B")

sim$run_up_to_time(sim$get_clock() + 30)

plot_tissue(sim, num_of_bins = 500)

sim$sample_cells("A2", bottom_left = c(480, 480), top_right = c(480 + bbox_width, 480 + bbox_width))
sim$sample_cells("B2", bottom_left = c(500, 500), top_right = c(500 + bbox_width, 500 + bbox_width))

plot_tissue(sim, num_of_bins = 500)

# sim$add_genotype(name = "C",
#                  epigenetic_rates = c("+-" = 0.05, "-+" = 0.1),
#                  growth_rates = c("+" = 0.7, "-" = 0.5),
#                  death_rates = c("+" = 0.05, "-" = 0.05))
# 
# cell <- sim$choose_cell_in("B")
# sim$mutate_progeny(cell, "C")
# 
# sim$update_rates("B-", death=1)
# sim$update_rates("B+", death=1)
# 
# sim$run_up_to_size("C+", 100)
# 
# plot_tissue(sim, num_of_bins = 500)

plot_sampled_cells(sim)

forest = sim$get_samples_forest()
cells = forest$get_nodes() %>% 
  filter(sample == 'A1', epistate == '-') %>% 
  pull(cell_id)

forest$get_coalescent_cells(cells)


set.seed(1234)

sim <- new(Simulation)
sim$update_tissue(250, 250)

sim$add_genotype(name = "A", growth_rates = 0.08, death_rates = 0.01)
sim$place_cell("A", 125, 125)
sim$run_up_to_size("A", 300)

plot_tissue(sim)

# Sampling with random box sampling
ncells = 20
n_w = n_h = 7

bbox = bbox_sampler(sim, "A", ncells, n_w, n_h)
sim$sample_cells("A1", bbox$p, bbox$q)

bbox = bbox_sampler(sim, "A", ncells, n_w, n_h)
sim$sample_cells("B1", bbox$p, bbox$q)

plot_tissue(sim)

# New subclone
sim$add_genotype(name = "B", growth_rates = 0.13, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size("B", 200)

plot_tissue(sim)

# Sampling with random box sampling
bbox = bbox_sampler(sim, "A", ncells, n_w, n_h)
sim$sample_cells("A2", bbox$p, bbox$q)

bbox = bbox_sampler(sim, "B", ncells, n_w, n_h)
sim$sample_cells("B2", bbox$p, bbox$q)

plot_tissue(sim)

# New subclone
sim$add_genotype(name = "C", growth_rates = 0.24, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("A"), "C")

sim$run_up_to_size("C", 500)
plot_tissue(sim)

# Sampling with random box sampling
bbox = bbox_sampler(sim, "B", ncells, n_w, n_h)
sim$sample_cells("A3", bbox$p, bbox$q)

bbox = bbox_sampler(sim, "C", ncells, n_w, n_h)
sim$sample_cells("B3", bbox$p, bbox$q)

plot_tissue(sim, num_of_bins = 100)

# Tree data
forest = sim$get_samples_forest()
forest_nodes = forest$get_nodes()

tree_plot = plot_sampled_cells(forest_nodes = forest_nodes) %>% 
  annotate_forest(forest)

tree_plot




plot_timeseries(sim)

plot_timeseries(sim) + ggplot2::scale_y_log10()

plot_muller(sim)

plot_muller(sim) + ggplot2::xlim(50, sim$get_clock()) + ggplot2::scale_y_log10()


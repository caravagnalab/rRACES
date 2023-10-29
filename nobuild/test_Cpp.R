devtools::load_all()

# Create a new simulation and set the simulation data directory
# to be "test". A 1000x1000-cells tissue is automatically
# added to the simulation
sim <- new(Simulation, "test")

# Get the simulation directory, i.e., "test"
sim$get_directory()

# Get the tissue size, i.e., c(1000,1000)
sim$get_tissue_size()

# Set the tissue name and size
sim$set_tissue("Liver", 1200, 900)

# Get the tissue name, i.e., "Liver"
sim$get_tissue_name()

# Get the tissue size, i.e., c(1200,900)
sim$get_tissue_size()

# Add two species, "A" and "B", having epigenetic states
sim$add_species(name = "A",
                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
                growth_rates = c("+" = 0.2, "-" = 0.08),
                death_rates = c("+" = 0.1, "-" = 0.01))

sim$add_species(name = "B",
                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
                growth_rates = c("+" = 0.25, "-" = 0.15),
                death_rates = c("+" = 0.05, "-" = 0.05))

# Program a timed transition from species "A" to
# species "B" at time 30
sim$add_timed_mutation("A", "B", 60)

# Get the simulation species with epigenetic state names
sim$get_species_names()

# Add one cell to the simulated tissue
sim$add_cell(genotype_name = "A+", x = 500, y = 500)

# Run the simulation up until the there are less than 100 
# cells in the species B-
sim$run_up_to_size("B-", 100)

# Get simulation time
sim$get_clock()

# Report the number of cells per species at the current simulation time, e.g.,
#   species epistate counts
# 1       A        -  30942
# 2       A        + 105715
# 3       B        -    100
# 4       B        +   1371
sim$get_counts()

# Run the simulation up to time 110
sim$run_up_to_time(110)

# Get simulation time, i.e., 110
sim$get_clock()

# Get the cells in the tissue at current simulation time
sim$get_cells()

# Get the cells in the tissue rectangular sample having [500,500]
# and [505,505] as lower and upper corners, respectively
sim$get_cells(c(500, 500), c(505, 505))

# Get the cells in the tissue having epigenetic state "-"
sim$get_cells(c("A", "B"), c("-"))

# Get the cells in the tissue having epigenetic state "-" and, at
# the same time, belonging to rectangular sample bounded by [500,500]
# and [505,505] as lower and upper corners, respectively
sim$get_cells(c(500, 500), c(505, 505), c("A", "B"), c("-"))

# plot_simulation_counts(sim)

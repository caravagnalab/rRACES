devtools::load_all()

# Create a new simulation and set the simulation name
# to be "test". A 1000x1000-cells tissue is automatically
# added to the simulation
sim <- new(Simulation, "test")

# Get the simulation name, i.e., "test"
sim$get_name()

# Get the tissue size, i.e., c(1000,1000)
sim$get_tissue_size()

# Set the tissue name and size
sim$update_tissue("Liver", 1200, 900)

# Get the tissue name, i.e., "Liver"
sim$get_tissue_name()

# Get the tissue size, i.e., c(1200,900)
sim$get_tissue_size()

# add genotype "A" together with their species
sim$add_genotype(name = "A",
                 epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
                 growth_rates = c("+" = 0.2, "-" = 0.08),
                 death_rates = c("+" = 0.1, "-" = 0.01))

# add genotype "B" together with their species
sim$add_genotype(name = "B",
                 epigenetic_rates = c("+-" = 0.05, "-+" = 0.05),
                 growth_rates = c("+" = 0.1, "-" = 0.3),
                 death_rates = c("+" = 0.05, "-" = 0.1))

# Schedule a mutation at time 60 from one of the species of
# genotype "A" to the species of "B" having the same
# epigenetic state
sim$schedule_genotype_mutation("A", "B", 60)

# Get the simulation species
sim$get_species()

# Add one cell to the simulated tissue
sim$add_cell(species = "A+", x = 500, y = 500)

# Run the simulation up until the there are less than 100
# cells in the species B-
sim$run_up_to_size("B-", 100)

# Get the number of fired event per species
sim$get_firings()

# Run the simulation up until number of deaths in B-
# is at least 100
sim$run_up_to_event("death", "B-", 100)

# Get the number of fired event per species. The row
# "death", "B", "-" reports 100 firings
sim$get_firings()

# Get simulation time
sim$get_clock()

# Report the number of cells per species at the current simulation time, e.g.,
#   species epistate counts
# 1       A        -  30942
# 2       A        + 105715
# 3       B        -    100
# 4       B        +   1371
# 5       C        -      0
# 6       C        +      0
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


# add genotype "C" together with their species
sim$add_genotype(name = "C",
                 epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
                 growth_rates = c("+" = 0.25, "-" = 0.15),
                 death_rates = c("+" = 0.05, "-" = 0.05))


# Plot the state of the simulation (number of cells)
plot_state(sim)

# Plot the cells in the tissue at current simulation time
plot_tissue(sim)

print(sim)
devtools::load_all()

# Create a new simulation and set the simulation data directory
# to be "output_dir". A 1000x1000-cells tissue is automatically
# added to the simulation
x <- create_simulation(
  name = "My Test rRACES",
  output_dir = "some_folder",
  tissue = "A fake liver",
  tissue_size = c(4000, 4000)
)

# Add two species, "A" and "B", having epigenetic states
x = add_species(
    x,
    name = "A",
    epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
    growth_rates = c("+" = 0.2, "-" = 0.08),
    death_rates = c("+" = 0.1, "-" = 0.01),
    ancestor = NULL,
    time_of_origin = 0
)

x = add_species(
  x,
  name = "B",
  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
  growth_rates = c("+" = 0.25, "-" = 0.15),
  death_rates = c("+" = 0.05, "-" = 0.05),
  ancestor = "A",
  time_of_origin = 10
)

# Program a timed transition from species "A" to
# species "B" at time 60
x$simulation$add_timed_mutation("A", "B", 60)

# Add one cell to the simulated tissue
x = set_initial_cell(
  x,
  species = "A",
  epistate = '+',
  position = c(500, 500)
)

# Run the simulation up to time 110
x = run(x, time = 110)

# Plot the state of the simulation (number of cells)
plot_state(x)


# Get the cells in the tissue at current simulation time
x$simulation$get_cells()

# Get the cells in the tissue rectangular sample having [500,500]
# and [505,505] as lower and upper corners, respectively
x$simulation$get_cells(c(500, 500), c(505, 505))

# Get the cells in the tissue having epigenetic state "-"
sim$get_cells(c("A", "B"), c("-"))

# Get the cells in the tissue having epigenetic state "-" and, at
# the same time, belonging to rectangular sample bounded by [500,500]
# and [505,505] as lower and upper corners, respectively
sim$get_cells(c(500, 500), c(505, 505), c("A", "B"), c("-"))

# Report the number of cells per species at the current simulation time
sim$get_counts()

# plot_simulation_counts(sim)

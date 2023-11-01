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
    death_rates = c("+" = 0.1, "-" = 0.01)
)

x = add_species(
  x,
  name = "B",
  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
  growth_rates = c("+" = 0.25, "-" = 0.15),
  death_rates = c("+" = 0.05, "-" = 0.05)
)

# Program a timed transition from species "A" to
# species "B" at time 60
x = add_genotype_evolution(x, "A", "B", 60)

# Add one cell to the simulated tissue
x = set_initial_cell(
  x,
  genotype = "A",
  epistate = "+",
  position = c(500, 500)
)

# Run the simulation up to some time
x = run(x, what = list(time = 80))

# Get the counts for species A+
Ap_count = x$simulation$get_counts() %>% 
  filter(genotype == "A", epistate == '+') %>% 
  pull(counts)

print(Ap_count)

# Run until we get 1000 more of A+
x = run(x, what = list(genotype = 'A', epistate = '+', size = Ap_count + 1000))

# Check that we have Ap_count + 1000
x$simulation$get_counts() %>% 
  filter(genotype == "A", epistate == '+') %>% 
  pull(counts)

# See how many firings we had
x$simulation$get_firings()

# Get the counts for species A+
Ap_switches = x$simulation$get_firings() %>% 
  filter(genotype == "A", epistate == '+', events == 'switch') %>% 
  pull(firings)

print(Ap_switches)

# Run until we get 200 more of A+ switching to A-
x = run(x, what = list(genotype = 'A', epistate = '+', event = 'switch', count = Ap_switches + 200))

# Check that we have Ap_switches + 200
x$simulation$get_firings() %>% 
  filter(genotype == "A", epistate == '+', events == 'switch') %>% 
  pull(firings)

x$simulation$get_counts() %>% 
  filter(genotype == "A", epistate == '+') %>% 
  pull(counts)

# Plot the state of the simulation (number of cells)
plot_state(x)

# Plot the cells in the tissue at current simulation time
plot_tissue(x)

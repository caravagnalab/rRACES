# useDynLib(rRACES, .registration=TRUE)
# import(methods, Rcpp)
# exportPattern("^[[:alpha:]]+")

devtools::load_all()

# create a new simulation and save simulation data in "output_dir"
sim <- new(Simulation, "output_dir")

# add two species
sim$add_species(name="A",
                epigenetic_rates=c("+-"=0.01,"-+"=0.01),
                growth_rates=c("+"=0.2,"-"=0.08),
                death_rates=c("+"=0.1,"-"=0.01))

sim$add_species(name="B",
                epigenetic_rates=c("+-"=0.01,"-+"=0.01),
                growth_rates=c("+"=0.1,"-"=0.3),
                death_rates=c("+"=0.05,"-"=0.1))

# Create a species object with the required parameters
# species <- list(
#   name = "B",
#   epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
#   growth_rates = c("+" = 0.1, "-" = 0.3),
#   death_rates = c("+" = 0.05, "-" = 0.1)
# )
# 
# # Call the R wrapper function with the simulation object and the prepared species object
# add_species(sim, species)


# get the simulation genotype names
sim$get_genotype_names()

# add one cell to the simulation
sim$add_cell(genotype_name="A+", x=500, y=500)

# simulate 60 time units
sim$run_up_to(60)

# plot_simulation_counts(sim)

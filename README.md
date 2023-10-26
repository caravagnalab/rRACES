# rRACES

rRACES is a R wrapper for [RACES](https://github.com/albertocasagrande/RACES).


## Compiling and Installing

```shell
git clone https://github.com/caravagnalab/rRACES.git
R CMD build rRACES
R CMD install rRACES_*.tar.gz
```

## Using `rRACES`

rRACES is at early stage of development and misses many features of 
RACES. At the moment, it supports:
- creating a spacial simulation
- adding species to the simulation
- adding cells to the spacial simulation
- running a simulation up to a specified simulated time

Further classes and methods will be available soon.

## A Small Example

```R
require(rRACES)

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

# get the simulation genotype names
sim$get_genotype_names()

# add one cell to the simulation
sim$add_cell(genotype_name="A+", x=500, y=500)

# simulate 60 time units
sim$run_up_to(60)
```

At the end of the computation, the directory `output_dir` will contain 
simulated data usable by RACES. 

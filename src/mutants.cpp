/*
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023-2024 Alberto Casagrande <alberto.casagrande@uniud.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <string>

#include <Rcpp.h>

#include "simulation.hpp"
#include "samples_forest.hpp"
#include "tissue_rectangle.hpp"

using namespace Rcpp;

namespace RE = Races::Mutants::Evolutions;
namespace RC = Races::Mutants;

RCPP_MODULE(Mutants){

//' @name TissueRectangle
//' @title A rectangle in the tissue
//' @field get_lower_corner Get the rectangle lower corner
//' @field get_upper_corner Get the rectangle upper corner
  class_<TissueRectangle>("TissueRectangle")
//' @name TissueRectangle$new
//' @title Build a new rectangle of tissue.
//' @examples
//' # build the rectangle [500,550]x[450,475]
//' rect <- new(TissueRectangle, c(500, 450), c(550, 475))
//'
//' rect
//'
//' # build the rectangle [500,550]x[450,475]
//' rect <- new(TissueRectangle, c(500, 450), 50, 25)
//'
//' rect
  .constructor<std::vector<uint16_t>, std::vector<uint16_t>>("Create a new rectangle")
  .constructor<std::vector<uint16_t>, RE::AxisSize, RE::AxisSize>("Create a new rectangle")

//' @name TissueRectangle$lower_corner
//' @title The lower corner of the tissue rectangle.
//' @examples
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the simulation death activation level
//' rect$lower_corner
  .property("lower_corner",&TissueRectangle::get_lower_corner, "The rectangle lower corner")

//' @name TissueRectangle$upper_corner
//' @title The lower corner of the tissue rectangle.
//' @examples
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the simulation death activation level
//' rect$upper_corner
  .property("upper_corner",&TissueRectangle::get_upper_corner, "The rectangle upper corner")
  .method("show",&TissueRectangle::show);

//' @name Simulation
//' @title Simulates the cell evolution on a tissue
//' @description The objects of this class can simulate the evolution
//'   of many cells belonging to different *species* on a tissue. Each
//'   cell can duplicate or die according to the rates that delineate
//'   the cell species.
//'
//'   `Simulation` supports epigenetic evolutions, and it lets users
//'   define species pairs that belong to the same mutant (even though,
//'   its genomic characterization is unknown) and differ because
//'   of their epigenetic state (i.e., either "+" or "-").
//'
//'   `Simulation` models epigenetic mutations and allows a cell in
//'   one of mutant species to generate a new cell belonging to
//'   the other species of the same mutant at a specified rate.
//'
//'   `Simulation` also allows users to schedule mutations from one
//'   mutant to a different mutant.
//' @field add_mutant Adds a mutant and its species \itemize{
//' \item \emph{Parameter:} \code{mutant} - The mutant name.
//' \item \emph{Parameter:} \code{epigenetic_rates} - The epigenetic rates of the mutant species (optional).
//' \item \emph{Parameter:} \code{growth_rates} - The duplication rates of the mutant species.
//' \item \emph{Parameter:} \code{death_rates} - The death rates of the mutant species.
//' }
//' @field choose_cell_in Chooses one cell in a mutant \itemize{
//' \item \emph{Parameter:} \code{mutant} - The mutant of the cell to choose.
//' \item \emph{Parameter:} \code{lower_corner} - The lower left corner of a rectangular selection (optional).
//' \item \emph{Parameter:} \code{upper_corner} - The upper right corner of a rectangular selection (optional).
//' \item \emph{Returns:} A list reporting "cell_id", "mutant", "epistate", "position_x",
//'   and "position_y" of the choosen cell.
//' }
//' @field death_activation_level The number of cells that activates cell death in a species.
//' @field duplicate_internal_cells Enable/disable duplication for internal cells.
//' @field get_added_cells Gets the cells manually added to the simulation \itemize{
//' \item \emph{Returns:} A data frame reporting "mutant", "epistate", "position_x",
//'   "position_y", and "time" for each cells manually added to
//'   the simulation.
//' }
//' @field search_sample Seach a rectangular sample having a minimum number of cells\itemize{
//' \item \emph{Parameter:} \code{mutant_name} - The mutant of the searched cells.
//' \item \emph{Parameter:} \code{num_of_cells} - The number of cells in the searched sample.
//' \item \emph{Parameter:} \code{width} - The width of the searched sample.
//' \item \emph{Parameter:} \code{height} - The height of the searched sample.
//' \item \emph{Returns:} If a rectangular sample satisfying the provided constraints can
//'   be found, the corresponding rectangle.
//' }
//' @field get_cell Gets one the tissue cells \itemize{
//' \item \emph{Parameter:} \code{x} - The position of the aimed cell on the x axis.
//' \item \emph{Parameter:} \code{y} - The position of the aimed cell on the y axis.
//' \item \emph{Returns:} A data frame reporting "cell_id", "mutant", "epistate", "position_x",
//'   and "position_y" of the aimed cell.
//' }
//' @field get_cells Gets the tissue cells \itemize{
//' \item \emph{Parameter:} \code{lower_corner} - The lower-left corner of the selection frame (optional).
//' \item \emph{Parameter:} \code{upper_corner} - The upper-right corner of the selection frame (optional).
//' \item \emph{Parameter:} \code{mutant_filter} - The vector of the to-be-selected mutant names (optional).
//' \item \emph{Parameter:} \code{epigenetic_filter} - The vector of the to-be-selected epigenetic states (optional).
//' \item \emph{Returns:} A data frame reporting "cell_id", "mutant", "epistate", "position_x",
//'   and "position_y" for each cells satisfying the provided filters and laying
//'   in the input frame.
//' }
//' @field get_clock Gets the simulated time \itemize{
//' \item \emph{Returns:} The time simulated by the simulation.
//' }
//' @field get_count_history Gets the history of the number of cells per species \itemize{
//' \item \emph{Returns:} A data frame reporting "mutant", "epistate", "counts",
//'   and "time" for each species and for each sampled time.
//' }
//' @field get_counts Counts the number of cells \itemize{
//' \item \emph{Returns:} A data frame reporting "mutant", "epistate", "counts" for each
//'   species in the simulation.
//' }
//' @field get_firing_history Gets the history of the number of fired events \itemize{
//' \item \emph{Returns:} A data frame reporting "event", "mutant", "epistate", "fired",
//'   and "time" for each event type, for each species, and for each sampled time.
//' }
//' @field get_firings Gets the number of fired events \itemize{
//' \item \emph{Returns:} A data frame reporting "event", "mutant", "epistate", and "fired"
//'   for each event type and for each species.
//' }
//' @field get_name Gets the simulation name \itemize{
//' \item \emph{Returns:} The simulation name, which corresponds to the name of the directory
//'   in which the simulation is saving its progresses.
//' }
//' @field get_lineage_graph Gets the simulation lineage graph\itemize{
//' \item \emph{Returns:} A data frame reporting "ancestor", "progeny", and "first_occurrence"
//'   of each species-to-species transition.
//' }
//' @field get_rates Gets the rates of a species\itemize{
//' \item \emph{Parameter:} \code{species} - The species whose rates are aimed.
//' \item \emph{Returns:} The list of the species names.
//' }
//' @field get_samples_forest Get the samples forest\itemize{
//' \item \emph{Returns:} The descendants forest having as leaves the sampled cells.
//' }
//' @field get_samples_info Retrieve information about the samples \itemize{
//' \item \emph{Returns:} A data frame containing, for each sample collected
//'   during the simulation, the columns "name", "time", "ymin",
//'   "xmin", "ymax", "xmax", and  "tumoral cells". "ymin",
//'   "xmin", "ymax", "xmax" report the boundaries of the sampled
//'   rectangular region, while "tumoral cells" is the number of
//'   tumoral cells in the sample.
//' }
//' @field get_species Gets the species \itemize{
//' \item \emph{Returns:} A data frame describing the registered species.
//' }
//' @field get_tissue_name Gets the tissue name \itemize{
//' \item \emph{Returns:} The name of the simulated tissue.
//' }
//' @field get_tissue_size Gets the size of the simulated tissue \itemize{
//' \item \emph{Returns:} The vector `c(x_size, y_size)` of the simulated tissue.
//' }
//' @field mutate_progeny Generate a mutated offspring \itemize{
//' \item \emph{Parameter:} \code{cell_position} - The position of the cell whose offspring will mutate.
//' \item \emph{Parameter:} \code{mutated_mutant} - The mutant of the mutated cell.
//' }
//' or
//' \itemize{
//' \item \emph{Parameter:} \code{x} - The position of the cell whose progeny will mutate on the x axis.
//' \item \emph{Parameter:} \code{y} - The position of the cell whose progeny will mutate on the y axis.
//' \item \emph{Parameter:} \code{mutated_mutant} - The mutant of the mutated cell.
//' }
//' @field place_cell Place one cell in the tissue \itemize{
//' \item \emph{Parameter:} \code{species} - The name of the new cell species.
//' \item \emph{Parameter:} \code{x} - The position on the x axis of the cell.
//' \item \emph{Parameter:} \code{y} - The position on the y axis of the cell.
//' }
//' @field schedule_mutation Schedules a mutant mutation \itemize{
//' \item \emph{Parameter:} \code{src} - The name of the mutant from which the mutation occurs.
//' \item \emph{Parameter:} \code{dest} - The name of the mutant to which the mutation leads.
//' \item \emph{Parameter:} \code{time} - The simulated time at which the mutation will occurs.
//' }
//' @field run_up_to_event Simulates cell evolution \itemize{
//' \item \emph{Parameter:} \code{event} - The considered event type, i.e., "growth", "death", or "switch".
//' \item \emph{Parameter:} \code{species} - The species whose event number is considered.
//' \item \emph{Parameter:} \code{num_of_events} - The threshold for the event number.
//' }
//' @field run_up_to_size Simulates cell evolution \itemize{
//' \item \emph{Parameter:} \code{species} - The species whose number of cells is considered.
//' \item \emph{Parameter:} \code{num_of_cells} - The threshold for the cell number.
//' }
//' @field run_up_to_time Simulates cell evolution \itemize{
//' \item \emph{Parameter:} \code{time} - The final simulation time.
//' }
//' @field sample_cells Sample a tissue rectangle region \itemize{
//' \item \emph{Parameter:} \code{name} - The sample name.
//' \item \emph{Parameter:} \code{lower_corner} - The bottom-left corner of the rectangle.
//' \item \emph{Parameter:} \code{upper_corner} - The top-right corner of the rectangle.
//' }
//' @field update_rates Updates the rates of a species\itemize{
//' \item \emph{Parameter:} \code{species} - The species whose rates must be updated.
//' \item \emph{Parameter:} \code{rates} - The list of the rates to be updated.
//' \item \emph{Returns:} The vector of the species names.
//' }
//' @field update_tissue Updates tissue name and size \itemize{
//' \item \emph{Parameter:} \code{name} - The new name of the tissue (optional).
//' \item \emph{Parameter:} \code{width} - The width of the new tissue.
//' \item \emph{Parameter:} \code{height} - The height of the new tissue.
//' }
//' @field var Builds a variable representing a simulation quantity \itemize{
//' \item \emph{Parameter:} \code{variable_description} - The description of
//'    the variable to be built.
//' \item \emph{Returns:} A variable representing the simulation quantity
//'   according to the parameter `variable_description`.
//' }
  class_<Simulation>("Simulation")

//' @name Simulation$new
//' @title Constructs a new Simulation
//' @param simulation_name The name of the simulation (optional).
//' @param seed The seed for the pseudo-random generator (optional).
//' @param save_snapshots A flag to save simulation snapshots on disk (optional,
//'   default `FALSE`).
//' @examples
//' # create a Simulation object storing binary dump in a temporary directory.
//' # The data are deleted from the disk as soon as the object is destroyed.
//' sim <- new(Simulation, "test")
//'
//' # add a new species, place a cell in the tissue, and let the simulation evolve.
//' sim$add_mutant(name = "A", growth_rate = 0.3, death_rate = 0.02)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_time(30)
//'
//' # no directory "test" has been created
//' "test" %in% list.files(".")
//'
//' # By using the optional parameter `save_snapshots`, we force the
//' # simulation to save its progresses in a local directory whose name
//' # is the name of the simulation, i.e., "test". This data will be
//' # preserved when the simulation object will be destroyed.
//' sim <- new(Simulation, "test", save_snapshots=TRUE)
//'
//' # as done above, we add a new species, place a cell in the tissue, and let the
//' # simulation evolve.
//' sim$add_mutant(name = "A", growth_rate = 0.3, death_rate = 0.02)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_time(30)
//'
//' # the directory "test" exists and contains a binary dump of
//' # sthe simulation.
//' "test" %in% list.files(".")
//'
//' # let us manually delete the "test" directory
//' unlink("test", recursive = TRUE)
//'
//' # we can also provide a random seed to the simulation...
//' sim <- new(Simulation, "test", 13)
//'
//' # ...or creating a simulation without providing any name. By default, the
//' # simulation name will have the following format `races_<date>_<hour>`.
//' sim <- new(Simulation, 13)
  .constructor("Create a simulation whose name is \"races_<year>_<hour><minute><second>\"")
  .constructor<SEXP>("Create a simulation")
  .constructor<SEXP, SEXP>("Crete a simulation")
  .constructor<std::string, int, bool>("Crete a simulation: the first parameter is the simulation name; "
                                       "the second one is the random seed; the third one is a Boolean flag "
                                       "to enable/disable simulation saves")

//' @name Simulation$place_cell
//' @title Place one cell in the tissue
//' @param species The name of the new cell species.
//' @param x The position on the x axis of the cell.
//' @param y The position on the y axis of the cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # add into the tissue a cell of species "A+" in position (500,500)
//' sim$place_cell("A+", 500, 500)
  .method("place_cell", &Simulation::place_cell, "Place a cell in the tissue")

//' @name Simulation$add_mutant
//' @title Adds a mutant and its species
//' @description This method adds a mutant and its species to the
//'   simulation. If the optional parameter `epigenetic_rate` is
//'   provided, then two new species having the same mutant and
//'   opposite epigenetic states are created. When, instead, the
//'   optional parameter `epigenetic_rate` is missing, this
//'   method creates only one species with no epigenetic states.
//' @param mutant The mutant name.
//' @param epigenetic_rates The epigenetic rates of the mutant species (optional).
//' @param growth_rates The duplication rates of the mutant species.
//' @param death_rates The death rates of the mutant species.
//' @examples
//' sim <- new(Simulation)
//'
//' # create the two species "A+" and "A-". They both have mutant "A".
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # create the species "C" its mutant is "C".
//' sim$add_mutant(name = "C", growth_rate = 0.2, death_rate = 0.1)
  .method("add_mutant", (void (Simulation::*)(const std::string&, const List&, const List&,
                                                const List&))(&Simulation::add_mutant),
          "Add a new species with epigenetic status")
  .method("add_mutant", (void (Simulation::*)(const std::string&, const double&,
                                                const double&))(&Simulation::add_mutant),
          "Add a new species")

//' @name Simulation$choose_cell_in
//' @title Chooses one cell in a mutant
//' @description This method chooses one of the cells whose mutant
//'   is `mutant`. Optionally, the lower and upper corners
//'   of a tissue rectangular selection can be provided
//'   to obtain one cell in the rectangle.
//' @param mutant The mutant of the cell to choose.
//' @param lower_corner The lower corner of the rectangular selection (optional).
//' @param upper_corner The upper corner of the rectangular selection (optional).
//' @return A list reporting `cell_id`, `mutant`, `epistate`, `position_x`,
//'   and `position_y` of the choosen cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.1, "-+" = 0.01),
//'                growth_rates = c("+" = 0.15, "-" = 0.3),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$death_activation_level <- 100
//' sim$schedule_mutation("A","B",20)
//' sim$run_up_to_size(species = "B-", num_of_cells = 50)
//'
//' # Randomly choose one cell in "B" in the tissue
//' sim$choose_cell_in(mutant = "B")
  .method("choose_cell_in", (List (Simulation::*)(const std::string&))(&Simulation::choose_cell_in),
          "Randomly choose one cell in a mutant")

  .method("choose_cell_in", (List (Simulation::*)(const std::string&,
                                                  const std::vector<RE::AxisPosition>&,
                                                  const std::vector<RE::AxisPosition>&))(&Simulation::choose_cell_in),
          "Randomly choose one cell having a specified mutant in a rectangular selection")

//' @name Simulation$schedule_mutation
//' @title Schedules a mutant mutation
//' @description This method schedules a mutant mutation that can occur
//'   from any of the species of the source mutant to the species of
//'   the destination mutant with a consistent epigenetic state.
//'   For the sake of example, if the mutation from "A" to "B" is
//'   scheduled, then we have three possible situations:
//'   1. The mutant "A" consists of the only species "A". Then,
//'      during one duplication of a cell of "A", one cell of "B"
//'      will arise.
//'   2. The mutant "A" consists of the species "A+" and "A-" and
//'      during one duplication of a cell of "A+", one cell of "B+"
//'      will arise.
//'   3. The mutant "A" consists of the species "A+" and "A-" and
//'      during one duplication of a cell of "A-", one cell of "B-"
//'      will arise.
//'   No other scenario can occur.
//' @param src The name of the mutant from which the mutation occurs.
//' @param dest The name of the mutant to which the mutation leads.
//' @param time The simulated time at which the mutation will occurs.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                growth_rates = c("+" = 0.3, "-" = 0.1),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # schedule an evolution from mutant "A" to mutant "B" at time 50
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
  .method("schedule_mutation", &Simulation::schedule_mutation,
          "Add a timed mutation between two different species")

//' @name Simulation$get_species
//' @title Gets the species
//' @return A data frame reporting `mutant`, `epistate`, `growth_rate`,
//'   `death_rate`, and `switch_rate` for each registered species.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_mutant("B", growth_rate = 0.15, death_rate = 0.05)
//'
//' # get the added species and their rates. In this case, "A"
//' # and "B"
//' sim$get_species()
  .method("get_species", &Simulation::get_species,
          "Get the species added to the simulation")

//' @name Simulation$get_samples_forest
//' @title Get the samples forest
//' @return The samples forest having as leaves the sampled cells
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' forest
  .method("get_samples_forest", &Simulation::get_samples_forest,
          "Get the descendants forest having as leaves the sampled cells")

//' @name Simulation$death_activation_level
//' @title The number of cells that activates cell death in a species.
//' @description This value is the minimum number of cells that
//'   enables cell death in a species. The cell of a species $S$ can die
//'   if and only if that $S$ has reached the death activation level at
//'   least once during the simulation.
//' @examples
//' sim <- new(Simulation)
//'
//' # get the simulation death activation level
//' sim$death_activation_level
//'
//' # set the death activation level to 50
//' sim$death_activation_level <- 50
  .property("death_activation_level", &Simulation::get_death_activation_level,
                                      &Simulation::set_death_activation_level,
            "The number of cells in a species that activates cell death" )

//' @name Simulation$duplicate_internal_cells
//' @title Enable/disable duplication for internal cells.
//' @description This Boolean flag enable/disable duplication of internal
//'   cells. When it is set to `FALSE`, the border-growth model
//'   is used. Otherwise, the homogeneous-growth model is applied.
//'   It is set to `FALSE` by default.
//' @examples
//' sim <- new(Simulation)
//'
//' # is the duplication of internal cells enabled? (by default, no)
//' sim$duplicate_internal_cells
//'
//' # enable homogeneous-growth model
//' sim$duplicate_internal_cells <- TRUE
//'
//' # now it should be set to `TRUE`
//' sim$duplicate_internal_cells
//'
//' # enable boder-growth model
//' sim$duplicate_internal_cells <- FALSE
  .property("duplicate_internal_cells", &Simulation::get_duplicate_internal_cells,
                                        &Simulation::set_duplicate_internal_cells,
            "Enable/disable duplication for internal cells" )

//' @name Simulation$get_clock
//' @title Gets the simulated time
//' @return The time simulated by the simulation.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the simulated time
//' sim$get_clock()
  .method("get_clock", &Simulation::get_clock, "Get the current simulation time")

//' @name Simulation$get_cell
//' @title Gets one of the tissue cells
//' @description This method collects some data of the aimed cell without altering
//'   the tissue.
//' @param x The position of the aimed cell on the x axis.
//' @param y The position of the aimed cell on the y axis.
//' @return A data frame reporting `cell_id`, `mutant`, `epistate`, `position_x`,
//'   and `position_y` of the aimed cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.02, "-" = 0.01))
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                growth_rates = c("+" = 0.3, "-" = 0.1),
//'                death_rates = c("+" = 0.02, "-" = 0.01))
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(40)
//'
//' # collect all the cells in the tissue
//' sim$get_cell(501, 502)
  .method("get_cell", (List (Simulation::*)(const RE::AxisPosition&,
                                            const RE::AxisPosition&) const)(&Simulation::get_cell),
          "Get one cell from the simulated tissue")

//' @name Simulation$get_cells
//' @title Gets the tissue cells
//' @description This method collects some data about the cells in the tissue
//'   without altering the tissue itself. The pairs of optional parameters
//'   `lower_corner` and `upper_corner` define a frame of the tissue in
//'   which the data are sampled. The optional parameters `mutant_filter`
//'   and `epigenetic_filter` filter the collected cell data according to
//'   the cell mutant and epigenetic state.
//' @param lower_corner The lower-left corner of the selection frame (optional).
//' @param upper_corner The upper-right corner of the selection frame (optional).
//' @param mutant_filter The vector of the to-be-selected mutant names (optional).
//' @param epigenetic_filter The vector of the to-be-selected epigenetic states (optional).
//' @return A data frame reporting `cell_id`, `mutant`, `epistate`, `position_x`,
//'   and `position_y` for each cells satisfying the provided filters and laying
//'   in the input frame.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                growth_rates = c("+" = 0.3, "-" = 0.1),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(30)
//'
//' # collect all the cells in the tissue
//' sim$get_cells()
//'
//' # get the cells in the frame [495,505]x[490,500]
//' sim$get_cells(lower_corner=c(495,490), upper_corner=c(505,500))
//'
//' # cells can be filtered by mutant name...
//' sim$get_cells(mutant_filter=c("A"),epigenetic_filter=c("+","-"))
//'
//' # ...or by epigenetic state
//' sim$get_cells(mutant_filter=c("A","B"),epigenetic_filter=c("-"))
//'
//' # cells can be filtered by frame, mutant, and epigenetic states
//' sim$get_cells(lower_corner=c(495,495), upper_corner=c(505,505),
//'               mutant_filter=c("A"),epigenetic_filter=c("+","-"))
  .method("get_cells", (List (Simulation::*)(const std::vector<RE::AxisPosition>&,
                                             const std::vector<RE::AxisPosition>&,
                                             const std::vector<std::string>&,
                                             const std::vector<std::string>&) const)(&Simulation::get_cells),
          "Get cells from the simulated tissue")
  .method("get_cells", (List (Simulation::*)(const SEXP&, const SEXP&) const)(&Simulation::get_cells),
          "Get cells from the simulated tissue")
  .method("get_cells", (List (Simulation::*)() const)(&Simulation::get_cells),
          "Get cells from the simulated tissue")

//' @name Simulation$get_name
//' @title Gets the simulation name
//' @return The simulation name, which corresponds to the name of the directory
//'   in which the simulation is saving its progresses.
//' @examples
//' sim <- new(Simulation)
//'
//' # Expecting "test"
//' sim$get_name()
  .method("get_name", &Simulation::get_name, "Get the simulation name")

//' @name Simulation$get_lineage_graph
//' @title Gets the simulation lineage graph
//' @description At the beginning of the computation only the species of the added
//'   cells are present in the tissue. As the simulation proceeds new species
//'   arise as a consequence of either mutant mutations or epigenetic
//'   switches. The *lineage graph* stores these species evolutions and it
//'   reports the first occurrence time of any species-to-species transition.
//'
//'   This method returns the lineage graph of the simulation.
//' @return A data frame reporting `ancestor`, `progeny`, and `first_occurrence` of
//'   each species-to-species transition.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                growth_rates = c("+" = 0.3, "-" = 0.1),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_mutation(src = "A", dst = "B", time = 20)
//' sim$run_up_to_time(50)
//'
//' sim$get_lineage_graph()
  .method("get_lineage_graph", &Simulation::get_lineage_graph,
          "Get the simulation lineage graph")

//' @name Simulation$get_tissue_name
//' @title Gets the tissue name
//' @return The name of the simulated tissue.
//' @examples
//' sim <- new(Simulation)
//' sim$update_tissue("Liver", 1200, 900)
//'
//' # get the tissue name, i.e., expecting "Liver"
//' sim$get_tissue_name()
  .method("get_tissue_name", &Simulation::get_tissue_name, "Get the simulation tissue name")

//' @name Simulation$get_tissue_size
//' @title Gets the size of the simulated tissue
//' @return The vector `c(x_size, y_size)` of the simulated tissue.
//' @examples
//' sim <- new(Simulation)
//' sim$update_tissue("Liver", 1200, 900)
//'
//' # get the tissue size, i.e., expecting c(1200,900)
//' sim$get_tissue_size()
  .method("get_tissue_size", &Simulation::get_tissue_size, "Get the simulation tissue size")

//' @name Simulation$get_added_cells
//' @title Gets the cells manually added to the simulation
//' @return A data frame reporting `mutant`, `epistate`, `position_x`,
//'   `position_y`, and `time` for each cells manually added to
//'   the simulation.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                growth_rates = c("+" = 0.3, "-" = 0.1),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_mutation(src = "A", dst = "B", time = 30)
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(50)
//'
//' # counts the number of cells per species
//' sim$get_added_cells()
  .method("get_added_cells", &Simulation::get_added_cells,
          "Get the cells manually added to the simulation")

//' @name Simulation$get_counts
//' @title Counts the number of cells
//' @return A data frame reporting `mutant`, `epistate`, `counts` for each
//'   species in the simulation.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_mutant("B", growth_rate = 0.15, death_rate = 0.05)
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_counts()
  .method("get_counts", &Simulation::get_counts, "Get the current number of cells per species")

//' @name Simulation$get_count_history
//' @title Gets the history of the number of cells per species
//' @description This method returns a data frame reporting the number of
//'   species cells in each sampled simulation time.
//' @return A data frame reporting `mutant`, `epistate`, `counts`,
//'   and `time` for each species, and for each sampled time.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_mutant("B", growth_rate = 0.15, death_rate = 0.05)
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A", 500, 500)
//' sim$history_delta <- 20
//' sim$run_up_to_time(70)
//'
//' # get the history of species counts
//' sim$get_count_history()
  .method("get_count_history", (List (Simulation::*)() const)&Simulation::get_count_history,
          "Get the number of simulated events per species along the computation")

//' @name Simulation$get_firings
//' @title Gets the number of fired events
//' @return A data frame reporting `event`, `mutant`, `epistate`, and `fired`
//'   for each event type, mutant, and epigenetic states.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the number of event fired per event and species
//' sim$get_firings()
  .method("get_firings", &Simulation::get_firings,
          "Get the current number of simulated events per species")

//' @name Simulation$get_firing_history
//' @title Gets the history of the number of fired events
//' @description This method returns a data frame reporting the number of
//'   events fired up to each sampled simulation time.
//' @return A data frame reporting `event`, `mutant`, `epistate`, `fired`,
//'   and `time` for each event type, for each species, and for each
//'   sampled time.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$history_delta <- 20
//' sim$run_up_to_time(70)
//'
//' # get the number of event fired per event and species
//' sim$get_firing_history()
  .method("get_firing_history", (List (Simulation::*)() const)&Simulation::get_firing_history,
          "Get the number of simulated events per species along the computation")

//' @name Simulation$get_rates
//' @title Get the rates of a species
//' @param species The species whose rates are aimed.
//' @return The list of the species rates.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Get the rates of "A-". In this case c("growth"=0.08, "death"=0.01, "switch"=0.02) is expected
//' sim$get_rates("A-")
  .method("get_rates", &Simulation::get_rates,
          "Get the rates of a species")

//' @name Simulation$get_samples_info
//' @title Retrieve information about the samples
//' @description This method retrieves information about
//'   the samples collected along the simulation.
//'   It returns a data frame reporting, for each
//'   sample, the name, the sampling time, the
//'   position, and the number of tumoural cells.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475),
//'                  upper_corner=c(500,550))
//'
//' # simulate 1 time unit more
//' sim$run_up_to_time(sim$get_clock()+1)
//'
//' # sample the region [500,520]x[525,550]
//' sim$sample_cells("S2", lower_corner=c(500,525),
//'                  upper_corner=c(520,550))
//'
//' # get information about all the collected
//' # samples, i.e, S1 and S2
//' sim$get_samples_info()
  .method("get_samples_info", &Simulation::get_samples_info,
          "Get some pieces of information about the collected samples")

//' @name Simulation$history_delta
//' @title The delta time between time series samples
//' @description This value is the maximum time between two successive
//'   time series data samples.
//' @examples
//' sim <- new(Simulation)
//'
//' # get the delta time between two time series samples (0 by default)
//' sim$history_delta
//'
//' # set the delta time between two time series samples
//' sim$death_activation_level <- 20
  .property("history_delta", &Simulation::get_history_delta,
                             &Simulation::set_history_delta,
            "The sampling delta for the get_*_history functions" )

//' @name Simulation$mutate_progeny
//' @title Generate a mutated progeny
//' @description This method simulates both the duplication of the cell in the
//'   specified position and the birth of one cells of a given
//'   mutant that preserves the epigenetic status of the original cell.
//'   The mutated cell will be located in the position of its parent.
//' @param cell_position The position of the cell whose offspring will mutate.
//' @param mutated_mutant The mutant of the mutated cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.01, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(30)
//'
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.1, "-+" = 0.01),
//'                growth_rates = c("+" = 0.15, "-" = 0.3),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # duplicate the cell in position (503, 492). One of
//' # its direct descendents will have mutant "B"
//' # sim$mutate_progeny(503, 492, "B")
//'
//' # the output of `choose_cell_in` and `get_cell` can also be used
//' # as input for `mutate_progeny`
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
  .method("mutate_progeny",  (void (Simulation::*)(const List&, const std::string&))
                                                  (&Simulation::mutate_progeny),
          "Duplicate a cell and mutate one of its children")
  .method("mutate_progeny",  (void (Simulation::*)(const RE::AxisPosition&,
                                                   const RE::AxisPosition&, const std::string&))
                                                  (&Simulation::mutate_progeny),
          "Duplicate a cell and mutate one of its children")

//' @name Simulation$run_up_to_time
//' @title Simulates cell evolution
//' @param time The final simulation time.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$place_cell("A", 500, 500)
//'
//' # simulate the tissue up to simulate timed 100
//' sim$run_up_to_time(40)
  .method("run_up_to_time", &Simulation::run_up_to_time,
          "Simulate the system up to the specified simulation time")

//' @name Simulation$run_up_to_event
//' @title Simulates cell evolution
//' @description This method simulates cell evolution until the number of events that
//'   have occurred to cells of a species reaches a specified threshold.
//' @param event The considered event, i.e., `growth`, `death`, or `switch`.
//' @param species The species whose event number is considered.
//' @param num_of_events The threshold for the event number.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//'
//' # simulate the cell evolution until the number of epigenetic events from
//' # the species "A+" is less than 100.
//' sim$run_up_to_event(event = "switch", species = "A+", num_of_events = 100)
  .method("run_up_to_event", &Simulation::run_up_to_event,
          "Simulate the system up to the specified number of events")

//' @name Simulation$run_up_to_size
//' @title Simulates cell evolution
//' @description This method simulates cell evolution until the number of cells in
//'   a species reaches a specified threshold.
//' @param species The species whose number of cells is considered.
//' @param num_of_cells The threshold for the cell number.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//'
//' # simulate the tissue until the species "A+" account for 100
//' # contemporary cells
//' sim$run_up_to_size(species = "A+", num_of_cells = 100)
  .method("run_up_to_size", &Simulation::run_up_to_size,
          "Simulate the system up to the specified number of cells in the species")

//' @name Simulation$run_until
//' @title Simulates cell evolution until a formula does not hold
//' @description This method simulates cell evolution until a formula does not
//'    hold.
//' @param formula The formula that will be satisfied at the end of the
//'    simulation.
//' @seealso `Simulation$var()`
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//'
//' # get the variable representing the simulation time
//' v_time <- sim$var("Time")
//'
//' # get the variable representing the cardinality of A+
//' va_p <- sim$var("A+")
//'
//' # get the variable representing the cardinality of A-
//' va_m <- sim$var("A-")
//'
//' # get the variable representing the number of epigenetic
//' # switches from A+
//' va_ps <- sim$var("A+.switches")
//'
//' # build a condition stating that the cardinality of A+ doubles
//' # that of A-
//' c1 <- va_p >= 2*va_m
//' 
//' # build a condition that holds when there are more than
//' # 100000 live cells of mutant A
//' c2 <- va_p + va_m > 1e5
//' 
//' # build a condition that holds when less than 4000 switched
//' # from A+ have occured
//' c3 <- va_ps < 4000
//' 
//' # build a condition that holds when 40 time unit have been
//' # simulated at least
//' c4 <- v_time >= 40
//' 
//' # build a condition that holds when c4 and at least one
//' # among c1, c2, and c3 hold
//' c5 <- c4 & (c1 | c2 | c3)
//' c5
//' 
//' # run the simulation while c5 does not hold
//' sim$run_until(c5)
//' 
//' sim
//' sim$get_clock()
  .method("run_until", &Simulation::run_until,
          "Simulate the system until a formula is *not* satisfied")

//' @name Simulation$sample_cells
//' @title Sample a tissue rectangle region.
//' @description This method removes a rectangular region from the simulated
//'   tissue and stores its cells in a sample that can subsequently
//'   retrieved to build a samples forest.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
  .method("sample_cells", &Simulation::sample_cells,
          "Sample a rectangular region of the tissue")

//' @name Simulation$update_rates
//' @title Update the rates of a species
//' @param species The species whose rates must be updated.
//' @param rates The list of rates to be updated.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Set the death and epigenetic switch rates of "A-" to 0
//' sim$update_rates("A-", c(switch=0, death=0))
  .method("update_rates", &Simulation::update_rates,
          "Update the rates of a species")

//' @name Simulation$update_tissue
//' @title Update tissue name and size
//' @param name The new name of the tissue (optional).
//' @param width The width of the new tissue.
//' @param height The height of the new tissue.
//' @examples
//' sim <- new(Simulation)
//'
//' # set the tissue size, but not the name
//' sim$update_tissue(1200, 900)
//'
//' # set the tissue size and its name
//' sim$update_tissue("Liver", 1200, 900)
  .method("update_tissue", (void (Simulation::*)(const std::string&, const RE::AxisSize&,
                                                 const RE::AxisSize&))(&Simulation::update_tissue),
          "Update tissue name and size")
  .method("update_tissue", (void (Simulation::*)(const RE::AxisSize&,
                                                 const RE::AxisSize&))(&Simulation::update_tissue),
          "Update tissue size")

//' @name Simulation$search_sample
//' @title Search a rectangular sample containing a minimum number of cells
//' @description This method searches a rectangular tissue sample containing
//'   the provided number of cells. The sizes of the sample are also
//'   provided a parameter of the method.
//'   The complexity of this method is O(|tissue rows|*|tissue cols|).
//' @param min_num_of_cells A named integer vector reporting the minimum number
//'   of cells per species or mutant.
//' @param num_of_cells The number of cells in the searched sample.
//' @param width The width of the searched sample.
//' @param height The height of the searched sample.
//' @return If a rectangular sample satisfying the provided constraints can
//'   be found, the corresponding rectangle.
//' @examples
//' sim <- new(Simulation)
//' sim$death_activation_level <- 50
//' sim$add_mutant(name = "A", growth_rate = 0.2, death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_size(species = "A", num_of_cells = 500)
//'
//' sim$add_mutant(name = "B", growth_rate = 0.3, death_rate = 0.01)
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
//' sim$run_up_to_size(species = "B", num_of_cells = 1000)
//'
//' # find a 10x10 sample containing 80 "B" cells and 10 "A" cells at least
//' sim$search_sample(c("A" = 10, "B" = 80),50,50)
  .method("search_sample", &Simulation::search_sample,
          "Search a rectangular sample containing a given number of cells")

//' @name Simulation$var
//' @title Builds a variable representing a simulation quantity
//' @description This method build a logic variable representing one of the
//'   simulation quantities among:
//'   -  cardinality of a species
//'   -  number of event among duplications, deaths, and epigenetic switches
//'   -  elapsed evolution time
//' @param variable_description The description of the variable to be built.
//'   When `variable_description` is the string `"Time"`, the elapsed
//'   simulation time variable is returned. If `variable_description` is
//'   set to a species name, then the variable representing the cardinality
//'   of the species is built. Finally, when the parameter is a species name
//'   followed by `.` and one among `duplications`, `deaths`, or `switches`,
//'   the variable representing the number of event of the specified type
//'   occurred since the computation beginning in the species.
//' @return A variable representing the simulation quantity according to
//'   the parameter `variable_description`.
//' @seealso `Simulation$run_until()`
//' @examples
//' # build a simulation and add two species to it
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # get the variable representing the simulation time
//' sim$var("Time")
//'
//' # get the variable representing the cardinality of A+
//' sim$var("A+")
//'
//' # get the variable representing the cardinality of A-
//' sim$var("A-")
//'
//' # get the variable representing the number of epigenetic
//' # switches from A+
//' sim$var("A+.switches")
//'
//' # get the variable representing the number of duplications
//' # in A+
//' sim$var("A+.duplications")
//'
//' # get the variable representing the number of deaths in A+
//' sim$var("A+.deaths")
  .method("var", &Simulation::get_var,
          "Get a variable representing a simulation quantity");

//' @name recover_simulation
//' @title Load a simulation
//' @param name The name of the simulation to be recovered
//' @examples
//' # create a simulation having name "recover_simulation_test" and
//' # save its snapshots in a local directory
//' sim <- new(Simulation, "recover_simulation_test",
//'            save_snapshots=TRUE)
//'
//' # add the species of "A"
//' sim$add_mutant("A",
//'                epigenetic_rates=c("+-" = 0.01, "-+"=0.01),
//'                growth_rates = c("+"=0.1, "-"=0.01),
//'                death_rates = c("+"=0.05, "-"=0.005))
//'
//' # place a cell in the tissue
//' sim$place_cell("A+", 500, 500)
//'
//' # simulate up to time 50
//' sim$run_up_to_time(50)
//'
//' # show the simulation
//' sim
//'
//' # remove the object sim from the environment
//' rm(list=c("sim"))
//'
//' # the object pointed by sim does not exist any more
//' exists("sim")
//'
//' # recover the simulation from the directory "recover_simulation_test"
//' sim <- recover_simulation("recover_simulation_test")
//'
//' sim
//'
//' # delete dump directory
//' unlink("recover_simulation_test", recursive = TRUE)
  function("recover_simulation", &Simulation::load,
           "Recover a simulation");

//' @name SamplesForest
//' @title The forest of the sampled cell ancestors.
//' @description Represents the forest of the ancestors of the
//'   cells sampled during the computation. The leaves of
//'   this forest are the sampled cells.
//' @field get_coalescent_cells Retrieve most recent common ancestors\itemize{
//' \item \emph{Parameter:} \code{cell_ids} - The list of the identifiers of the
//'   cells whose most recent common ancestors are aimed (optional).
//' \item \emph{Return:} A data frame representing, for each of the identified
//'   cells, the identified (column "cell_id"), whenever the
//'   node is not a root, the ancestor identifier (column
//'   "ancestor"), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column "sample"), the mutant
//'   (column "mutant"), the epistate (column "epistate"),
//'   and the birth time (column "birth_time").
//' }
//' @field get_nodes Get the forest nodes \itemize{
//' \item \emph{Return:} A data frame representing, for each node
//'   in the forest, the identified (column "id"),
//'   whenever the node is not a root, the ancestor
//'   identifier (column "ancestor"), whenever the node
//'   was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column "sample"), the mutant (column
//'   "mutant"), the epistate (column "epistate"),
//'   and the birth time (column "birth_time").
//' }
//' @field get_samples_info Retrieve information about the samples \itemize{
//' \item \emph{Returns:} A data frame containing, for each sample collected
//'   during the simulation, the columns "name", "time", "ymin",
//'   "xmin", "ymax", "xmax", and  "tumoral cells". "ymin",
//'   "xmin", "ymax", "xmax" report the boundaries of the sampled
//'   rectangular region, while "tumoral cells" is the number of
//'   tumoral cells in the sample.
//' }
//' @field get_species_info Gets the species data\itemize{
//' \item \emph{Returns:} A data frame reporting "mutant" and "epistate"
//'   for each registered species.
//' }
//' @field get_sticks Compute the forest sticks \itemize{
//' \item \emph{Returns:} The list of the forest sticks. Each stick is represented as
//'   the list of cell identifiers labelling the nodes in the stick
//'   from the higher to the deeper in the forest.
//' }
//' @field get_subforest_for Build a subforest using as leaves some of the original samples \itemize{
//' \item \emph{Parameter:} \code{sample_names} - The names of the samples whose cells will be used
//'   as leaves of the new forest.
//' \item \emph{Returns:} A samples forest built on the samples mentioned in `sample_names`.
//' }
//' @field save Save a samples forest in a file \itemize{
//' \item \emph{Parameter:} \code{filename} - The path of the file in which the samples
//'   forest must be saved.
//' }
  class_<SamplesForest>("SamplesForest")

//' @name SamplesForest$get_nodes
//' @title Get the nodes of the forest
//' @return A data frame representing, for each node
//'   in the forest, the identified (column "cell_id"),
//'   whenever the node is not a root, the ancestor
//'   identifier (column "ancestor"), whenever the
//'   node was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column "sample"), the mutant (column
//'   "mutant"), the epistate (column "epistate"),
//'   and the birth time (column "birth_time").
//' @examples
//' # create a simulation
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' forest$get_nodes()
    .method("get_nodes", (List (SamplesForest::*)() const)(&SamplesForest::get_nodes),
            "Get the nodes of the forest")

//' @name SamplesForest$get_coalescent_cells
//' @title Retrieve most recent common ancestors
//' @description This method retrieves the most recent common ancestors
//'   of a set of cells. If the optional parameter `cell_ids` is
//'   used, this method find the most recent common ancestors of
//'   the cells having an identifier among those in `cell_ids`.
//'   If, otherwise, the optional parameter is not used, this
//'   method find the most recent common ancestors of the forest
//'   leaves.
//' @param cell_ids The list of the identifiers of the cells whose
//'   most recent common ancestors are aimed (optional).
//' @return A data frame representing, for each of the identified
//'   cells, the identified (column "cell_id"), whenever the
//'   node is not a root, the ancestor identifier (column
//'   "ancestor"), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column "sample"), the mutant
//'   (column "mutant"), the epistate (column "epistate"),
//'   and the birth time (column "birth_time").
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' forest$get_coalescent_cells()
    .method("get_coalescent_cells",
            (List (SamplesForest::*)(const std::list<Races::Mutants::CellId>&) const)
                (&SamplesForest::get_coalescent_cells),
            "Get the most recent common ancestor of some cells")
    .method("get_coalescent_cells",
            (List (SamplesForest::*)() const)(&SamplesForest::get_coalescent_cells),
            "Get the most recent common ancestor of all the forest trees")

//' @name SamplesForest$get_subforest_for
//' @title Build a subforest using as leaves some of the original samples
//' @param sample_names The names of the samples whose cells will be used
//'   as leaves of the new forest
//' @return A samples forest built on the samples mentioned in `sample_names`
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A",
//'                growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' sim$run_up_to_size(species = "A", num_of_cells = 60000)
//'
//' # sample again the same region
//' sim$sample_cells("S2", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' forest$get_subforest_for("S2")
    .method("get_subforest_for", &SamplesForest::get_subforest_for,
            "Get the sub-forest for some of the original samples")

//' @name SamplesForest$get_samples_info
//' @title Retrieve information about the samples
//' @description This method retrieves information about
//'   the samples whose cells were used as leaves
//'   of the samples forest.
//' @return A data frame reporting, for each sample, the
//'   name, the sampling time, the position, and
//'   the number of tumoural cells.
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A", growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' # get information about the sampled whose cells
//' # are the forest leaves, i.e, S1 and S2
//' forest$get_samples_info()
    .method("get_samples_info", &SamplesForest::get_samples_info,
            "Get some pieces of information about the samples")

//' @name SamplesForest$get_species_info
//' @title Gets the species
//' @return A data frame reporting `mutant` and `epistate`
//'   for each registered species.
    .method("get_species_info", &SamplesForest::get_species_info,
            "Get the recorded species")

//' @name SamplesForest$get_sticks
//' @title Compute the forest sticks
//' @description A _crucial node_ of a forest is a root of the forest, a node
//'   whose parent belongs to a different mutant, or the most recent
//'   common ancestor of two crucial nodes.
//'
//'   A _stick_ is a path of the forest in which the only crucial
//'   nodes are the first and the last one.
//'
//'   This method return the list of the forest sticks. Each stick is
//'   represented by the sequence of cell identifiers labelling the
//'   nodes in the stick.
//' @return The list of the forest sticks. Each stick is represented as
//'   the list of cell identifiers labelling the nodes in the stick
//'   from the higher to the deeper in the forest.
//' @seealso [PhylogeneticForest$get_sticks()]
//' @examples
//' sim <- new(Simulation)
//' sim$add_mutant(name = "A", growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 5000)
//'
//' sim$add_mutant(name = "B", growth_rate = 0.3, death_rate = 0.01)
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
//' sim$run_up_to_size(species = "B", num_of_cells = 1000)
//'
//' # search for a 33x33 region containing 50 cells in A and
//' # 50 cells in B at least and sample it
//' region <- sim$search_sample(c(A = 50, B = 50), 33, 33)
//' sim$sample_cells("S1", region$lower_corner, region$upper_corner)
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' # search for the forest sticks
//' forest$get_sticks()
    .method("get_sticks", (std::list<std::list<Races::Mutants::CellId>> (SamplesForest::*)() const)(&SamplesForest::get_sticks),
            "Get the forest sticks")

//' @name SamplesForest$save
//' @title Save a samples forest in a file
//' @param filename The path of the file in which the samples
//'   forest must be saved.
    .method("save", &SamplesForest::save,
            "Save a samples forest")

    .method("show", &SamplesForest::show,
            "Describe the SamplesForest");

//' @name load_samples_forest
//' @title Load a samples forest from a file
//' @param filename The path of the file from which the samples
//'   forest must be load.
//' @return The load samples forest
  function("load_samples_forest", &SamplesForest::load,
           "Recover a samples forest");
}
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

namespace RE = RACES::Mutants::Evolutions;
namespace RC = RACES::Mutants;

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
//' @title The tissue rectangle lower corner
//' @description This is the vector of the minima among all the rectangle dimensions.
//' @examples
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the simulation death activation level
//' rect$lower_corner
  .property("lower_corner",&TissueRectangle::get_lower_corner, "The rectangle lower corner")

//' @name TissueRectangle$upper_corner
//' @title The tissue rectangle upper corner
//' @description This is the vector of the maxima among all the rectangle dimensions.
//' @examples
//' rect <- new(TissueRectangle, c(500, 500), c(550, 550))
//'
//' # get the simulation death activation level
//' rect$upper_corner
  .property("upper_corner",&TissueRectangle::get_upper_corner, "The rectangle upper corner")
  .method("show",&TissueRectangle::show);

//' @name SpatialSimulation
//' @title Simulating the cell evolution on a tissue
//' @description This class simulates the cell evolution on a tissue.
//' @details The objects of this class can simulate the evolution
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
//' \item \emph{Returns:} A list reporting "`cell_id`", "`mutant`", "`epistate`", "`position_x`",
//'   and "`position_y`" of the choosen cell.
//' }
//' @field death_activation_level The number of cells that activates cell death in a species.
//' @field border_growth_model Switch between homogeneous and border driven growth models.
//' @field get_added_cells Gets the cells manually added to the simulation \itemize{
//' \item \emph{Returns:} A dataframe reporting "`mutant`", "`epistate`", "`position_x`",
//'   "`position_y`", and "`time`" for each cells manually added to
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
//' \item \emph{Returns:} A dataframe reporting "`cell_id`", "`mutant`", "`epistate`", "`position_x`",
//'   and "`position_y`" of the aimed cell.
//' }
//' @field get_cells Gets the tissue cells \itemize{
//' \item \emph{Parameter:} \code{lower_corner} - The lower-left corner of the selection frame (optional).
//' \item \emph{Parameter:} \code{upper_corner} - The upper-right corner of the selection frame (optional).
//' \item \emph{Parameter:} \code{mutant_filter} - The vector of the to-be-selected mutant names (optional).
//' \item \emph{Parameter:} \code{epigenetic_filter} - The vector of the to-be-selected epigenetic states (optional).
//' \item \emph{Returns:} A dataframe reporting "`cell_id`", "`mutant`", "`epistate`", "`position_x`",
//'   and "`position_y`" for each cells satisfying the provided filters and laying
//'   in the input frame.
//' }
//' @field get_clock Gets the simulated time \itemize{
//' \item \emph{Returns:} The time simulated by the simulation.
//' }
//' @field get_count_history Gets the history of the number of cells per species \itemize{
//' \item \emph{Returns:} A dataframe reporting "`mutant`", "`epistate`", "`counts`",
//'   and "`time`" for each species and for each sampled time.
//' }
//' @field get_counts Counts the number of cells per species\itemize{
//' \item \emph{Returns:} A dataframe reporting "`mutant`", "`epistate`", "`counts`", and
//'   "`overall`" for each species in the simulation.
//' }
//' @field get_firing_history Gets the history of the number of fired events \itemize{
//' \item \emph{Returns:} A dataframe reporting "event", "`mutant`", "`epistate`", "`fired`",
//'   and "`time`" for each event type, for each species, and for each sampled time.
//' }
//' @field get_firings Gets the number of fired events \itemize{
//' \item \emph{Returns:} A dataframe reporting "event", "`mutant`", "`epistate`", and "`fired`"
//'   for each event type and for each species.
//' }
//' @field get_name Gets the simulation name \itemize{
//' \item \emph{Returns:} The simulation name, which corresponds to the name of the directory
//'   in which the simulation is saving its progresses.
//' }
//' @field get_lineage_graph Gets the simulation lineage graph\itemize{
//' \item \emph{Returns:} A dataframe reporting "`ancestor`", "`progeny`", and "`first_occurrence`"
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
//' \item \emph{Returns:} A dataframe containing, for each sample collected
//'   during the simulation, the columns "`name`", "`time`", "`ymin`",
//'   "`xmin`", "`ymax`", "`xmax`", "`tumour_cells`", and
//'   "`tumour_cells_in_bbox`". The columns "`ymin`", "`xmin`", "`ymax`",
//'   "`xmax`" report the boundaries of the sample bounding box, while
//'   "`tumour_cells`" and "`tumour_cells_in_bbox`" are the number of tumour
//'   cells in the sample and in the bounding box, respectively.
//' }
//' @field get_species Gets the species \itemize{
//' \item \emph{Returns:} A dataframe describing the registered species.
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
//' \item \emph{Parameter:} \code{quiet} - A Boolean flag to avoid progress bar (optional).
//' }
//' @field run_up_to_size Simulates cell evolution \itemize{
//' \item \emph{Parameter:} \code{species} - The species whose number of cells is considered.
//' \item \emph{Parameter:} \code{num_of_cells} - The threshold for the cell number.
//' \item \emph{Parameter:} \code{quiet} - A Boolean flag to avoid progress bar (optional).
//' }
//' @field run_up_to_time Simulates cell evolution \itemize{
//' \item \emph{Parameter:} \code{time} - The final simulation time.
//' \item \emph{Parameter:} \code{quiet} - A Boolean flag to avoid progress bar (optional).
//' }
//' @field run_until Simulating cell evolution \itemize{
//' \item \emph{Parameter:} \code{formula} - The formula that will be satisfied at the
//'    end of the simulation.
//' \item \emph{Parameter:} \code{quiet} - A Boolean flag to avoid progress bar (optional).
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
//' @field var Builds a variable representing a simulation quantity \itemize{
//' \item \emph{Parameter:} \code{variable_description} - The description of
//'    the variable to be built.
//' \item \emph{Returns:} A variable representing the simulation quantity
//'   according to the parameter `variable_description`.
//' }
  class_<SpatialSimulation>("SpatialSimulation")

//' @name SpatialSimulation$place_cell
//' @title Placing one cell in the tissue
//' @description This method places a cell in the tissue.
//' @param species The name of the new cell species.
//' @param x The position on the x axis of the cell.
//' @param y The position on the y axis of the cell.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a tissue simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # add into the tissue a cell of species "A+" in position (500,500)
//' sim$place_cell("A+", 500, 500)
  .method("place_cell", &SpatialSimulation::place_cell,
          "Place a cell in the tissue")

//' @name SpatialSimulation$add_mutant
//' @title Adding a mutant and its species
//' @description This method adds a mutant and its species to the
//'   simulation.
//' @details If the optional parameter `epigenetic_rate` is
//'   provided, then two new species having the same mutant and
//'   opposite epigenetic states are created. When, instead, the
//'   optional parameter `epigenetic_rate` is missing, this
//'   method creates only one species with no epigenetic states.
//' @param mutant The mutant name.
//' @param epigenetic_rates The epigenetic rates of the mutant species (optional).
//' @param growth_rates The duplication rates of the mutant species.
//' @param death_rates The death rates of the mutant species.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//'
//' # create the two species "A+" and "A-". They both have mutant "A".
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # create the species "C" its mutant is "C".
//' sim$add_mutant(name = "C", growth_rate = 0.2, death_rate = 0.1)
  .method("add_mutant", (void (SpatialSimulation::*)(const std::string&, const List&, const List&,
                                                const List&))(&SpatialSimulation::add_mutant),
          "Add a new species with epigenetic status")
  .method("add_mutant", (void (SpatialSimulation::*)(const std::string&, const double&,
                                                const double&))(&SpatialSimulation::add_mutant),
          "Add a new species")

//' @name SpatialSimulation$choose_cell_in
//' @title Picking one cell in a mutant
//' @description This method chooses one cell among those of a mutant.
//' @details It randomly chooses one of the cells of a mutant.
//'   Optionally, the lower and upper corners of a tissue rectangular
//'   selection can be provided to obtain one cell in the rectangle.
//' @param mutant The mutant of the cell to choose.
//' @param lower_corner The lower corner of the rectangular selection (optional).
//' @param upper_corner The upper corner of the rectangular selection (optional).
//' @return A list reporting `cell_id`, `mutant`, `epistate`, `position_x`,
//'   and `position_y` of the choosen cell.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_mutant(name = "B",
//'                epigenetic_rates = c("+-" = 0.1, "-+" = 0.01),
//'                growth_rates = c("+" = 0.15, "-" = 0.4),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$death_activation_level <- 100
//' sim$schedule_mutation("A","B",20)
//' sim$run_up_to_size(species = "B-", num_of_cells = 50)
//'
//' # Randomly choose one cell in "B" in the tissue
//' sim$choose_cell_in(mutant = "B")
  .method("choose_cell_in", (List (SpatialSimulation::*)(const std::string&))(&SpatialSimulation::choose_cell_in),
          "Randomly choose one cell in a mutant")
  .method("choose_cell_in", (List (SpatialSimulation::*)(const std::string&,
                                                  const std::vector<RE::AxisPosition>&,
                                                  const std::vector<RE::AxisPosition>&))(&SpatialSimulation::choose_cell_in),
          "Randomly choose one cell having a specified mutant in a rectangular selection")

//' @name SpatialSimulation$schedule_mutation
//' @title Scheduling a mutant mutation
//' @description This method schedules a mutant mutation
//' @details The mutation can occur from any of the species of
//'   the source mutant to the species of the destination mutant
//'   with a consistent epigenetic state.
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
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("schedule_mutation", &SpatialSimulation::schedule_mutation,
          "Add a timed mutation between two different species")

//' @name SpatialSimulation$get_species
//' @title Getting the species
//' @description This method returns the simulated species.
//' @return A dataframe reporting `mutant`, `epistate`, `growth_rate`,
//'   `death_rate`, and `switch_rate` for each registered species.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_mutant("B", growth_rate = 0.15, death_rate = 0.05)
//'
//' # get the added species and their rates. In this case, "A"
//' # and "B"
//' sim$get_species()
  .method("get_species", &SpatialSimulation::get_species,
          "Get the species added to the simulation")

//' @name SpatialSimulation$get_samples_forest
//' @title Getting the samples forest
//' @description This method returns the samples forest.
//' @return The samples forest having as leaves the sampled cells
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("get_samples_forest", &SpatialSimulation::get_samples_forest,
          "Get the descendants forest having as leaves the sampled cells")

//' @name SpatialSimulation$death_activation_level
//' @title Death activation level
//' @description The number of cells that activates cell death in a species.
//' @details This value is the minimum number of cells that
//'   enables cell death in a species. The cell of a species $S$ can die
//'   if and only if that $S$ has reached the death activation level at
//'   least once during the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//'
//' # get the simulation death activation level
//' sim$death_activation_level
//'
//' # set the death activation level to 50
//' sim$death_activation_level <- 50
  .property("death_activation_level", &SpatialSimulation::get_death_activation_level,
                                      &SpatialSimulation::set_death_activation_level,
            "The number of cells in a species that activates cell death" )

//' @name SpatialSimulation$border_growth_model
//' @title Internal cells duplication
//' @description This property switches between homogeneous and border driven
//'   growth models.
//' @details This Boolean flag switches between homogeneous and border driven
//'   growth models. When it is set to `TRUE`, the border-growth model
//'   is used. Otherwise, the homogeneous-growth model is applied.
//'   It is set to `TRUE` by default.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//'
//' # is the simulation using the border driven growth model
//' # (default: TRUE)
//' sim$border_growth_model
//'
//' # switch to the homogeneous-growth model
//' sim$border_growth_model <- FALSE
//'
//' # now it is set to `FALSE`
//' sim$border_growth_model
//'
//' # switch back to the border-growth model
//' sim$border_growth_model <- FALSE
  .property("border_growth_model", &SpatialSimulation::is_border_growth_model,
                                   &SpatialSimulation::set_border_growth_model,
            "Switch between homogeneous and border driven growth models." )

//' @name SpatialSimulation$get_clock
//' @title Getting the simulated time
//' @description This method returns the current simulation time.
//' @return The time simulated by the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the simulated time
//' sim$get_clock()
  .method("get_clock", &SpatialSimulation::get_clock, "Get the current simulation time")

//' @name SpatialSimulation$get_cell
//' @title Getting one of the tissue cells
//' @description This method collects some data of the aimed cell without altering
//'   the tissue.
//' @param x The position of the aimed cell on the x axis.
//' @param y The position of the aimed cell on the y axis.
//' @return A dataframe reporting `cell_id`, `mutant`, `epistate`, `position_x`,
//'   and `position_y` of the aimed cell.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("get_cell", (List (SpatialSimulation::*)(const RE::AxisPosition&,
                                            const RE::AxisPosition&) const)(&SpatialSimulation::get_cell),
          "Get one cell from the simulated tissue")

//' @name SpatialSimulation$get_cells
//' @title Getting the tissue cells
//' @description This method returns information about tumour tissue cells
//' @details It collects some data about the cells in the tissue
//'   without altering the tissue itself. The pairs of optional parameters
//'   `lower_corner` and `upper_corner` define a frame of the tissue in
//'   which the data are sampled. The optional parameters `mutant_filter`
//'   and `epigenetic_filter` filter the collected cell data according to
//'   the cell mutant and epigenetic state.
//' @param lower_corner The lower-left corner of the selection frame
//'   (optional).
//' @param upper_corner The upper-right corner of the selection frame
//'   (optional).
//' @param mutant_filter The vector of the to-be-selected mutant names
//'   (optional).
//' @param epigenetic_filter The vector of the to-be-selected epigenetic states
//'   (optional).
//' @return A dataframe reporting "`cell_id`", "`mutant`", "`epistate`",
//'   "`position_x`", and "`position_y`" for each cells satisfying the provided
//'   filters and laying in the input frame.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("get_cells",
          (List (SpatialSimulation::*)(const std::vector<RE::AxisPosition>&,
                                const std::vector<RE::AxisPosition>&,
                                const std::vector<std::string>&,
                                const std::vector<std::string>&) const)
                    (&SpatialSimulation::get_cells),
          "Get cells from the simulated tissue")
  .method("get_cells",
          (List (SpatialSimulation::*)(const SEXP&, const SEXP&) const)
                    (&SpatialSimulation::get_cells),
          "Get cells from the simulated tissue")
  .method("get_cells",
          (List (SpatialSimulation::*)(const std::string&) const)(
            &SpatialSimulation::get_cells),
          "Get cells from what was the simulated tissue status before a sampling")
  .method("get_cells", (List (SpatialSimulation::*)() const)(&SpatialSimulation::get_cells),
          "Get cells from the simulated tissue")

//' @name SpatialSimulation$get_name
//' @title Getting the simulation name
//' @description This method returns the simulation name
//' @return The simulation name, which corresponds to the name of the directory
//'   in which the simulation is saving its progresses.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//'
//' # Expecting "test"
//' sim$get_name()
  .method("get_name", &SpatialSimulation::get_name, "Get the simulation name")

//' @name SpatialSimulation$get_lineage_graph
//' @title Getting the simulation lineage graph
//' @description This method returns the lineage graph of the simulation.
//' @details At the beginning of the computation only the species of the
//'   added cells are present in the tissue. As the simulation proceeds new
//'   species arise as a consequence of either mutant mutations or epigenetic
//'   switches. The *lineage graph* stores these species evolutions and it
//'   reports the first occurrence time of any species-to-species transition.
//' @return A dataframe reporting `ancestor`, `progeny`, and `first_occurrence`
//'   of each species-to-species transition.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("get_lineage_graph", &SpatialSimulation::get_lineage_graph,
          "Get the simulation lineage graph")

//' @name SpatialSimulation$get_tissue_size
//' @title Getting the simulated tissue size
//' @description This method returns the size of the simulated tissue.
//' @return The vector `c(x_size, y_size)` of the simulated tissue.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation having size 1200x900
//' sim <- SpatialSimulation(width=1200, height=900)
//'
//' # get the tissue size, i.e., expecting c(1200,900)
//' sim$get_tissue_size()
//' @seealso `SpatialSimulation()`
  .method("get_tissue_size", &SpatialSimulation::get_tissue_size, "Get the simulation tissue size")

//' @name SpatialSimulation$get_added_cells
//' @title Getting the cells manually added to the simulation
//' @description This method returns the cells manually added to
//'   the simulation.
//' @return A dataframe reporting `mutant`, `epistate`, `position_x`,
//'   `position_y`, and `time` for each cells manually added to
//'   the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("get_added_cells", &SpatialSimulation::get_added_cells,
          "Get the cells manually added to the simulation")

//' @name SpatialSimulation$get_counts
//' @title Counting the cell number
//' @description This method returns the current number of cells per
//'   species and that since the simulation began.
//' @return A dataframe reporting `mutant`, `epistate`, `counts`, and
//'   `overall` for each species in the simulation.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_mutant("B", growth_rate = 0.15, death_rate = 0.05)
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_counts()
  .method("get_counts", &SpatialSimulation::get_counts,
          "Get the current number of cells and that "
          "throughout the entire simulation")

//' @name SpatialSimulation$get_count_history
//' @title Getting the history of the number of cells per species
//' @description This method returns a dataframe reporting the number of
//'   species cells in each sampled simulation time.
//' @return A dataframe reporting `mutant`, `epistate`, `counts`,
//'   and `time` for each species, and for each sampled time.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_mutant("B", growth_rate = 0.15, death_rate = 0.05)
//' sim$schedule_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A", 500, 500)
//' sim$history_delta <- 20
//' sim$run_up_to_time(70)
//'
//' # get the history of species counts
//' sim$get_count_history()
  .method("get_count_history", (List (SpatialSimulation::*)() const)&SpatialSimulation::get_count_history,
          "Get the number of simulated events per species along the computation")

//' @name SpatialSimulation$get_firings
//' @title Getting the number of fired events
//' @description This method returns a dataframe reporting the current
//'   number of simulated events per species.
//' @return A dataframe reporting `event`, `mutant`, `epistate`, and `fired`
//'   for each event type, mutant, and epigenetic states.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the number of event fired per event and species
//' sim$get_firings()
  .method("get_firings", &SpatialSimulation::get_firings,
          "Get the current number of simulated events per species")

//' @name SpatialSimulation$get_firing_history
//' @title Getting the fired event history
//' @description This method returns a dataframe reporting the number of
//'   events fired up to each sampled simulation time.
//' @return A dataframe reporting `event`, `mutant`, `epistate`, `fired`,
//'   and `time` for each event type, for each species, and for each
//'   sampled time.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("get_firing_history",
          (List (SpatialSimulation::*)() const)&SpatialSimulation::get_firing_history,
          "Get the number of simulated events per species along the computation")

//' @name SpatialSimulation$get_rates
//' @title Getting the species rates
//' @description This method return the rates of a species.
//' @param species The species whose rates are aimed.
//' @return The list of the species rates.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Get the rates of "A-". In this case c("growth"=0.08, "death"=0.01,
//' # "switch"=0.02) is expected
//' sim$get_rates("A-")
  .method("get_rates", &SpatialSimulation::get_rates,
          "Get the rates of a species")

//' @name SpatialSimulation$get_rates_update_history
//' @title Retrieving the rates update history
//' @description This method retrieves the simulation rates
//'   update history.
//' @return A dataframe containing the event rates updates. The dataframe
//'   contains the columns "`time`", "`mutant`", "`epistate`", "`event`",
//'   and "`rate`". Each row reports an update in the rate of an event
//'   in a species.
//' @seealso [SpatialSimulation$update_rates()], [SpatialSimulation$get_rates()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' sim$place_cell("A+", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A-", num_of_cells = 5000)
//'
//' # Set the death and epigenetic switch rates of "A-" to 0
//' sim$update_rates("A-", c(switch=0, death=0))
//'
//' sim$run_up_to_size(species = "A+", num_of_cells = 5000)
//'
//' # Set the death rate of "A+" to 0.5
//' sim$update_rates("A+", c(death=0.5))
//'
//' sim$run_up_to_time(sim$get_clock()+1)
//'
//' # get the rates update history
//' sim$get_rates_update_history()
  .method("get_rates_update_history", &SpatialSimulation::get_rates_update_history,
          "Get the rates update history")

//' @name SpatialSimulation$get_samples_info
//' @title Retrieving sample information
//' @description This method retrieves information about
//'   the samples collected along the simulation.
//' @return A dataframe containing, for each sample collected
//'   during the simulation, the columns "`name`", "`time`", "`id`",
//'   "`ymin`", "`xmin`", "`ymax`", "`xmax`", "`tumour_cells`", and
//'   "`tumour_cells_in_bbox`". The columns "`ymin`", "`xmin`", "`ymax`",
//'   "`xmax`" report the boundaries of the sample bounding box, while
//'   "`tumour_cells`" and "`tumour_cells_in_bbox`" are the number of tumour
//'   cells in the sample and in the bounding box, respectively.
//' @seealso [SpatialSimulation$sample_cells()],
//'   [SamplesForest$get_samples_info()],
//'   [PhylogeneticForest$get_samples_info()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
//' # simulate 1 time unit more
//' sim$run_up_to_time(sim$get_clock()+1)
//'
//' # sample the region [500,520]x[525,550]
//' sim$sample_cells("S2", lower_corner=c(500,525), upper_corner=c(520,550))
//'
//' # get information about all the collected
//' # samples, i.e, S1 and S2
//' sim$get_samples_info()
  .method("get_samples_info",
          (List (SpatialSimulation::*)() const)&SpatialSimulation::get_samples_info,
          "Get some pieces of information about the collected samples")

//' @name SpatialSimulation$history_delta
//' @title The delta time between time series samples
//' @description This value is the maximum time between two successive
//'   time series data samples.
//' @seealso `Simulation`
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//'
//' # get the delta time between two time series samples (0 by default)
//' sim$history_delta
//'
//' # set the delta time between two time series samples
//' sim$death_activation_level <- 20
  .property("history_delta", &SpatialSimulation::get_history_delta,
                             &SpatialSimulation::set_history_delta,
            "The sampling delta for the get_*_history functions" )

//' @name SpatialSimulation$mutate_progeny
//' @title Generating a mutated progeny
//' @description This method generates a mutated progeny.
//' @details It simulates both the duplication of the cell in the
//'   specified position and the birth of one cells of a given
//'   mutant that preserves the epigenetic status of the original cell.
//'   The mutated cell will be located in the position of its parent.
//' @param cell_position The position of the cell whose offspring will mutate.
//' @param mutated_mutant The mutant of the mutated cell.
//' @seealso `Simulation`, [SpatialSimulation$choose_cell_in()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("mutate_progeny",  (void (SpatialSimulation::*)(const List&, const std::string&))
                                                  (&SpatialSimulation::mutate_progeny),
          "Duplicate a cell and mutate one of its children")
  .method("mutate_progeny",  (void (SpatialSimulation::*)(const RE::AxisPosition&,
                                                   const RE::AxisPosition&, const std::string&))
                                                  (&SpatialSimulation::mutate_progeny),
          "Duplicate a cell and mutate one of its children")

//' @name SpatialSimulation$run_up_to_time
//' @title Simulates cell evolution
//' @param time The final simulation time.
//' @param quiet An optional  Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @seealso `Simulation`, [SpatialSimulation$run_up_to_event()],
//'    [SpatialSimulation$run_up_to_size()], [SpatialSimulation$run_until()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$place_cell("A", 500, 500)
//'
//' # simulate the tissue up to simulate timed 100
//' sim$run_up_to_time(40)
  .method("run_up_to_time",
          (void (SpatialSimulation::*)(const RACES::Time&, const bool))
          &SpatialSimulation::run_up_to_time,
          "Simulate the system up to the specified simulation time")
  .method("run_up_to_time",
          (void (SpatialSimulation::*)(const RACES::Time&))
          &SpatialSimulation::run_up_to_time,
          "Simulate the system up to the specified simulation time")

//' @name SpatialSimulation$run_up_to_event
//' @title Simulating cell evolution
//' @description This method simulates cell evolution until the number of events that
//'   have occurred to cells of a species reaches a specified threshold.
//' @param event The considered event, i.e., `growth`, `death`, or `switch`.
//' @param species The species whose event number is considered.
//' @param num_of_events The threshold for the event number.
//' @param quiet An optional  Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @seealso `Simulation`, [SpatialSimulation$run_up_to_time()],
//'    [SpatialSimulation$run_up_to_size()], [SpatialSimulation$run_until()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//'
//' # simulate the cell evolution until the number of epigenetic events from
//' # the species "A+" is less than 100.
//' sim$run_up_to_event(event = "switch", species = "A+", num_of_events = 100)
  .method("run_up_to_event",
          (void (SpatialSimulation::*)(const std::string&, const std::string&,
                                       const size_t&, const bool))
          &SpatialSimulation::run_up_to_event,
          "Simulate the system up to the specified number of events")
  .method("run_up_to_event",
          (void (SpatialSimulation::*)(const std::string&, const std::string&,
                                       const size_t&))
          &SpatialSimulation::run_up_to_event,
          "Simulate the system up to the specified number of events")

//' @name SpatialSimulation$run_up_to_size
//' @title Simulating cell evolution
//' @description This method simulates cell evolution until the number of cells in
//'   a species reaches a specified threshold.
//' @param species The species whose number of cells is considered.
//' @param num_of_cells The threshold for the cell number.
//' @param quiet An optional  Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @seealso `Simulation`, [SpatialSimulation$run_up_to_time()],
//'    [SpatialSimulation$run_up_to_event()], [SpatialSimulation$run_until()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//'
//' # simulate the tissue until the species "A+" account for 100
//' # contemporary cells
//' sim$run_up_to_size(species = "A+", num_of_cells = 100)
  .method("run_up_to_size",
          (void (SpatialSimulation::*)(const std::string&, const size_t&,
                                       const bool))
          &SpatialSimulation::run_up_to_size,
          "Simulate the system up to the specified number of cells in the species")
  .method("run_up_to_size",
          (void (SpatialSimulation::*)(const std::string&, const size_t&))
          &SpatialSimulation::run_up_to_size,
          "Simulate the system up to the specified number of cells in the species")

//' @name SpatialSimulation$run_until
//' @title Simulating cell evolution
//' @description This method simulates cell evolution until a formula does not
//'    hold.
//' @param formula The formula that will be satisfied at the end of the
//'    simulation.
//' @param quiet An optional  Boolean flag to avoid the progress bar
//'   (default: FALSE).
//' @seealso `Simulation`, [SpatialSimulation$var()], [SpatialSimulation$run_up_to_time()],
//'    [SpatialSimulation$run_up_to_event()], [SpatialSimulation$run_up_to_size()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("run_until",
          (void (SpatialSimulation::*)(const Logics::Formula&,
                                       const bool))
          &SpatialSimulation::run_until,
          "Simulate the system until a formula is *not* satisfied")
  .method("run_until",
          (void (SpatialSimulation::*)(const Logics::Formula&))
          &SpatialSimulation::run_until,
          "Simulate the system until a formula is *not* satisfied")

//' @name SpatialSimulation$sample_cells
//' @title Sampling a set of cells
//' @description This method samples a set of tumour cells.
//' @details It removes the cells from the simulated tissue and
//'   stores them in a sample that can be subsequently
//'   retrieved to build a samples forest.
//' @param sample_name The name of the sample.
//' @param lower_corner The lower corner of the sample bounding box (optional
//'   in pair with `upper_corner`).
//' @param upper_corner The upper corner of the sample bounding box (optional
//'   in pair with `lower_corner`).
//' @param num_of_cells The maximum number of tumour cells to collect
//'   (optional).
//' @seealso [SpatialSimulation$get_samples_info()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # randomly sample 50 tumour cells from the tissue
//' sim$sample_cells("S1", num_of_cells=50)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S2", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # randomly sample 50 tumour cells from the tissue region [500,550]x[500,550]
//' sim$sample_cells("S3", lower_corner=c(500,500), upper_corner=c(550,550),
//'                  num_of_cells=50)
//'
//' sim$get_samples_info()
  .method("sample_cells",
          (void (SpatialSimulation::*)(const std::string&,
                                const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                                const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                                const size_t& num_of_cells) const)(&SpatialSimulation::sample_cells),
          "Sample a rectangular region of the tissue")
  .method("sample_cells",
          (void (SpatialSimulation::*)(const std::string&,
                                const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                                const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner) const)(&SpatialSimulation::sample_cells),
          "Sample a rectangular region of the tissue")
  .method("sample_cells",
          (void (SpatialSimulation::*)(const std::string&,
                                const size_t& num_of_cells) const)(&SpatialSimulation::sample_cells),
          "Sample a rectangular region of the tissue")

//' @name SpatialSimulation$update_rates
//' @title Updating species rates
//' @description This method updates the rates of a species.
//' @param species The species whose rates must be updated.
//' @param rates The list of rates to be updated.
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A",
//'                epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                growth_rates = c("+" = 0.2, "-" = 0.08),
//'                death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Set the death and epigenetic switch rates of "A-" to 0
//' sim$update_rates("A-", c(switch=0, death=0))
  .method("update_rates", &SpatialSimulation::update_rates,
          "Update the rates of a species")

//' @name SpatialSimulation$search_sample
//' @title Searching a rectangular tissue sample
//' @description This method searches a rectangular tissue sample.
//' @details The aimed sample mush satisfy the specified number of cells.
//'   The sizes of the samples are also provided a parameter of the method.
//'   The complexity of this method is
//'   \eqn{O(|\textrm{tissue width}|*|\textrm{tissue height}|)}.
//' @param min_num_of_cells A named integer vector reporting the minimum number
//'   of cells per species or mutant.
//' @param num_of_cells The number of cells in the searched sample.
//' @param width The width of the searched sample.
//' @param height The height of the searched sample.
//' @return If a rectangular sample satisfying the provided constraints can
//'   be found, the corresponding rectangle.
//' @seealso `Simulation`
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$death_activation_level <- 50
//' sim$add_mutant(name = "A", growth_rate = 0.2, death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_size(species = "A", num_of_cells = 500)
//'
//' sim$add_mutant(name = "B", growth_rate = 0.3, death_rate = 0.01)
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
//' sim$run_up_to_size(species = "B", num_of_cells = 1000)
//'
//' # find a 50x50 sample containing 80 "B" cells and 10 "A" cells at least
//' sim$search_sample(c("A" = 10, "B" = 80), 50, 50)
  .method("search_sample", &SpatialSimulation::search_sample,
          "Search a rectangular sample containing a given number of cells")

//' @name SpatialSimulation$search_samples
//' @title Searching rectangular tissue samples
//' @description This method searches a set of rectangular tissue samples.
//' @details The aimed samples mush satisfy the specified number of cells.
//'   The sizes of the samples are also provided a parameter of the method.
//'   This method takes asymptotic time
//'   \eqn{\Theta(|\textrm{tissue width}|*|\textrm{tissue height}|)}.
//' @param min_num_of_cells A named integer vector reporting the minimum number
//'   of cells per species or mutant.
//' @param num_of_cells The number of cells in the searched sample.
//' @param width The width of the searched sample.
//' @param height The height of the searched sample.
//' @param n_samples The number of searched samples.
//' @param seed The seed of the random generator the select the samples
//'     among those satisfying the constraints (optional).
//' @return A vector of `n_samples` rectangular tissue samples that
//'     satisfy the aimed constraints.
//' @seealso `Simulation`
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$death_activation_level <- 50
//' sim$add_mutant(name = "A", growth_rate = 0.2, death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_size(species = "A", num_of_cells = 50)
//'
//' sim$add_mutant(name = "B", growth_rate = 0.3, death_rate = 0.01)
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
//' sim$run_up_to_size(species = "B", num_of_cells = 40000)
//'
//' plot <- plot_tissue(sim, num_of_bins = 1000)
//'
//' # find 3 50x50 samples containing 80 "B" cells and 100 "A" cells
//' # at least
//' bboxes <- sim$search_samples(c("A" = 100, "B" = 80), 50, 50,
//'                              n_samples=3)
//' bboxes
//'
//' # plot the bbox of the found samples
//' for (bbox in bboxes) {
//'   plot <- plot +
//'     ggplot2::geom_rect(xmin = bbox$lower_corner[1],
//'                        xmax = bbox$upper_corner[1],
//'                        ymin = bbox$lower_corner[2],
//'                        ymax = bbox$upper_corner[2],
//'                        fill = NA, color = "black")
//' }
//'
//' plot
  .method("search_samples", (std::vector<TissueRectangle> (SpatialSimulation::*)(const Rcpp::IntegerVector&,
                                                                          const uint16_t&, const uint16_t&,
                                                                          const size_t) const)(&SpatialSimulation::search_samples),
          "Search rectangular samples containing a given number of cells")
  .method("search_samples", (std::vector<TissueRectangle> (SpatialSimulation::*)(const Rcpp::IntegerVector&,
                                                                          const uint16_t&, const uint16_t&,
                                                                          const size_t, int) const)(&SpatialSimulation::search_samples),
          "Search rectangular samples containing a given number of cells")

//' @name SpatialSimulation$var
//' @title Building a simulation status variable
//' @description This method builds a simulation status variable.
//' @details This method builds a logic variable representing one of the
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
//' @seealso `Simulation`, [SpatialSimulation$run_until()]
//' @examples
//' # build a simulation and add two species to it
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
  .method("var", &SpatialSimulation::get_var,
          "Get a variable representing a simulation quantity");

//' @name recover_simulation
//' @title Loading a simulation
//' @description This method loads a simulation from the disk.
//' @param name The name of the simulation to be recovered.
//' @examples
//' # create a simulation having name "recover_simulation_test" and
//' # save its snapshots in a local directory
//' sim <- SpatialSimulation("recover_simulation_test",
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
  function("recover_simulation", &SpatialSimulation::load,
           "Recover a simulation");

//' @name SpatialSimulation
//' @title Building a new simulation
//' @description This method builds a new simulation.
//' @param name The name of the simulation (default:
//'     "`races_<year>_<hour><minute><second>`").
//' @param width The width of the simulated tissue (default: 1000).
//' @param height The height of the simulated tissue (default: 1000).
//' @param save_snapshots A flag to save simulation snapshots on disk
//'   (default: `FALSE`).
//' @param seed The seed for the pseudo-random generator (optional).
//' @examples
//' # create a SpatialSimulation object storing binary dump in a temporary
//' # directory. The data are deleted from the disk as soon as the object
//' # is destroyed.
//' sim <- SpatialSimulation("test")
//'
//' # add a new species, place a cell in the tissue, and let the simulation
//' # evolve.
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
//' sim <- SpatialSimulation("test", save_snapshots=TRUE)
//'
//' # as done above, we add a new species, place a cell in the tissue,
//' # and let the simulation evolve.
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
//' sim <- SpatialSimulation("test", seed=13)
//'
//' # ...or creating a simulation without providing any name. By default, the
//' # simulation name will have the following format `races_<date>_<hour>`.
//' sim <- SpatialSimulation(seed=13)
//'
//' # finally we can also specify the size of the simulated space
//' # by using the optional parameters `width` and `height`
//' sim <- SpatialSimulation(width=1200, height=900)
  function("SpatialSimulation", &SpatialSimulation::build_simulation,
           List::create(_["name"] = R_NilValue, _["width"]=1000, _["height"]=1000,
                        _["save_snapshots"] = false, _["seed"] = R_NilValue),
            "Create a spatial simulation");

//' @name SamplesForest
//' @title The sample cell ancestor forest
//' @description This class represents the forest of the ancestors of the
//'   cells sampled during the computation. The leaves of
//'   this forest are the sampled cells.
//' @field get_coalescent_cells Retrieve most recent common ancestors\itemize{
//' \item \emph{Parameter:} \code{cell_ids} - The list of the identifiers of the
//'   cells whose most recent common ancestors are aimed (optional).
//' \item \emph{Return:} A dataframe representing, for each of the identified
//'   cells, the identified (column "`cell_id`"), whenever the
//'   node is not a root, the ancestor identifier (column
//'   "`ancestor`"), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column "`sample`"), the mutant
//'   (column "`mutant`"), the epistate (column "`epistate`"),
//'   and the birth time (column "`birth_time`").
//' }
//' @field get_nodes Get the forest nodes \itemize{
//' \item \emph{Return:} A dataframe representing, for each node
//'   in the forest, the identified (column "`cell_id`"),
//'   whenever the node is not a root, the ancestor
//'   identifier (column "`ancestor`"), whenever the node
//'   was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column "`sample`"), the mutant (column
//'   "`mutant`"), the epistate (column "`epistate`"),
//'   and the birth time (column "`birth_time`").
//' }
//' @field get_samples_info Retrieve information about the samples \itemize{
//' \item \emph{Returns:} A dataframe containing, for each sample collected
//'   during the simulation, the columns "`name`", "`time`", "`ymin`",
//'   "`xmin`", "`ymax`", "`xmax`", "`tumour_cells`", and
//'   "`tumour_cells_in_bbox`". The columns "`ymin`", "`xmin`", "`ymax`",
//'   "`xmax`" report the boundaries of the sample bounding box, while
//'   "`tumour_cells`" and "`tumour_cells_in_bbox`" are the number of tumour
//'   cells in the sample and in the bounding box, respectively.
//' }
//' @field get_species_info Gets the species data\itemize{
//' \item \emph{Returns:} A dataframe reporting "`mutant`" and "`epistate`"
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
//' @title Getting forest nodes
//' @description This method builds a dataframe containing forest nodes.
//' @return A dataframe representing, for each node
//'   in the forest, the identified (column "`cell_id`"),
//'   whenever the node is not a root, the ancestor
//'   identifier (column "`ancestor`"), whenever the
//'   node was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column "`sample`"), the mutant (column
//'   "`mutant`"), the epistate (column "`epistate`"),
//'   and the birth time (column "`birth_time`").
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
//' nodes <- forest$get_nodes()
//'
//' head(nodes, 5)
    .method("get_nodes", (List (SamplesForest::*)() const)(&SamplesForest::get_nodes),
            "Get the nodes of the forest")

//' @name SamplesForest$get_coalescent_cells
//' @title Retrieving the most recent common ancestors
//' @description This method retrieves the most recent common ancestors
//'   of a set of cells.
//' @details If the optional parameter `cell_ids` is
//'   used, this method find the most recent common ancestors of
//'   the cells having an identifier among those in `cell_ids`.
//'   If, otherwise, the optional parameter is not used, this
//'   method find the most recent common ancestors of the forest
//'   leaves.
//' @param cell_ids The list of the identifiers of the cells whose
//'   most recent common ancestors are aimed (optional).
//' @return A dataframe representing, for each of the identified
//'   cells, the identified (column "`cell_id`"), whenever the
//'   node is not a root, the ancestor identifier (column
//'   "`ancestor`"), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column "`sample`"), the mutant
//'   (column "`mutant`"), the epistate (column "`epistate`"),
//'   and the birth time (column "`birth_time`").
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
            (List (SamplesForest::*)(const std::list<RACES::Mutants::CellId>&) const)
                (&SamplesForest::get_coalescent_cells),
            "Get the most recent common ancestor of some cells")
    .method("get_coalescent_cells",
            (List (SamplesForest::*)() const)(&SamplesForest::get_coalescent_cells),
            "Get the most recent common ancestor of all the forest trees")

//' @name SamplesForest$get_subforest_for
//' @title Building subforests
//' @description This method builds a subforest using as leaves some of the original
//'   samples.
//' @param sample_names The names of the samples whose cells will be used
//'   as leaves of the new forest
//' @return A samples forest built on the samples mentioned in `sample_names`
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
//' @title Retrieving samples' information
//' @description This method retrieves information about
//'   the samples whose cells were used as leaves
//'   of the samples forest.
//' @return A dataframe containing, for each sample collected
//'   during the simulation, the columns "`name`", "`time`", "`id`",
//'   "`ymin`", "`xmin`", "`ymax`", "`xmax`", "`tumour_cells`", and
//'   "`tumour_cells_in_bbox`". The columns "`ymin`", "`xmin`", "`ymax`",
//'   "`xmax`" report the boundaries of the sample bounding box, while
//'   "`tumour_cells`" and "`tumour_cells_in_bbox`" are the number of tumour
//'   cells in the sample and in the bounding box, respectively.
//' @seealso [PhylogeneticForest$get_samples_info()] for usage examples,
//'   [SpatialSimulation$sample_cells()], [SpatialSimulation$get_samples_info()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
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
//' @title Getting forest species
//' @description This method builds a dataframe containing information
//'      about the simulated species.
//' @return A dataframe reporting `mutant` and `epistate`
//'   for each registered species.
    .method("get_species_info", &SamplesForest::get_species_info,
            "Get the recorded species")

//' @name SamplesForest$get_sticks
//' @title Computing the forest sticks
//' @description This method computes the forest sticks.
//' @details A _crucial node_ of a forest is a root of the forest, a node
//'   whose parent belongs to a different species, or the most recent common
//'   ancestor of two crucial nodes.
//'
//'   A _stick_ is a path of the forest in which the only crucial nodes are
//'   the first and the last one.
//'
//'   This method return the list of the forest sticks. Each stick is
//'   represented by the sequence of cell identifiers labelling the nodes in
//'   the stick.
//' @param birth_threshold The maximum birth time for the cells associated
//'   to the returned sticks (optional).
//' @return The list of the forest sticks whose associated cells have
//'   birth time smaller than or equal to `birth_threshold`. Each stick is
//'   represented as the list of cell identifiers labelling the nodes in the
//'   stick from the higher to the deeper in the forest.
//' @seealso [PhylogeneticForest$get_sticks()]
//' @examples
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- SpatialSimulation()
//' sim$add_mutant(name = "A", growth_rate = 0.2,
//'                death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 100)
//'
//' sim$get_clock()
//'
//' sim$add_mutant(name = "B", growth_rate = 0.3, death_rate = 0.01)
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
//' sim$run_up_to_size(species = "B", num_of_cells = 120)
//'
//' sim$get_clock()
//'
//' sim$add_mutant(name = "C", growth_rate = 0.3, death_rate = 0.01)
//' sim$mutate_progeny(sim$choose_cell_in("B"), "C")
//' sim$run_up_to_size(species = "C", num_of_cells = 200)
//'
//' # search for a 33x33 region containing 50 cells in A and
//' # 50 cells in B at least and sample it
//' region <- sim$search_sample(c(A = 50, B = 50), 33, 33)
//' sim$sample_cells("S1", region$lower_corner, region$upper_corner)
//'
//' # search for a 33x33 region containing 50 cells in B and
//' # 50 cells in C at least and sample it
//' region <- sim$search_sample(c(B = 50, C = 50), 33, 33)
//' sim$sample_cells("S2", region$lower_corner, region$upper_corner)
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' # search for the forest sticks
//' forest$get_sticks()
//'
//' # search for the forest sticks whose corresponding cells have
//' # birth times 40 time units at most
//' forest$get_sticks(40)
    .method("get_sticks", (std::list<std::list<RACES::Mutants::CellId>> (SamplesForest::*)(const double) const)(&SamplesForest::get_sticks),
            "Get the forest sticks")
    .method("get_sticks", (std::list<std::list<RACES::Mutants::CellId>> (SamplesForest::*)() const)(&SamplesForest::get_sticks),
            "Get the forest sticks")

//' @name SamplesForest$save
//' @title Saving samples forests
//' @description This method saves a samples forest in a file.
//' @param filename The path of the file in which the samples
//'   forest must be saved.
    .method("save", &SamplesForest::save,
            "Save a samples forest")

    .method("show", &SamplesForest::show,
            "Describe the SamplesForest");

//' @name load_samples_forest
//' @title Loading samples forests
//' @description This method loads a samples forest in a file.
//' @param filename The path of the file from which the samples
//'   forest must be load.
//' @return The load samples forest
  function("load_samples_forest", &SamplesForest::load,
           "Recover a samples forest");
}
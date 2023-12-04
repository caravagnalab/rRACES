/*
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023 Alberto Casagrande <alberto.casagrande@uniud.it>
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

namespace RS = Races::Drivers::Simulation;
namespace RD = Races::Drivers;

RCPP_EXPOSED_CLASS(TissueRectangle)
RCPP_EXPOSED_CLASS(Simulation)
RCPP_EXPOSED_CLASS(SamplesForest)
RCPP_MODULE(Drivers){
  class_<TissueRectangle>("TissueRectangle")
  .constructor<std::vector<uint16_t>, std::vector<uint16_t>>("Create a new rectangle")
  .constructor<std::vector<uint16_t>, RS::AxisSize, RS::AxisSize>("Create a new rectangle")
  .property("lower_corner",&TissueRectangle::get_lower_corner, "The rectangle lower corner")
  .property("upper_corner",&TissueRectangle::get_upper_corner, "The rectangle upper corner")
  .method("show",&TissueRectangle::show);

  class_<Simulation>("Simulation")
  .constructor("Create a simulation whose name is \"races_<year>_<hour><minute><second>\"")
  .constructor<SEXP>("Create a simulation")
  .constructor<SEXP, SEXP>("Crete a simulation")
  .constructor<std::string, int, bool>("Crete a simulation: the first parameter is the simulation name; "
                                       "the second one is the random seed; the third one is a Boolean flag "
                                       "to enable/disable simulation saves")

  // place_cell
  .method("place_cell", &Simulation::place_cell, "Place a cell in the tissue")

  // add_genotype
  .method("add_genotype", (void (Simulation::*)(const std::string&, const List&, const List&,
                                                const List&))(&Simulation::add_genotype),
          "Add a new species with epigenetic status")
  .method("add_genotype", (void (Simulation::*)(const std::string&, const double&,
                                                const double&))(&Simulation::add_genotype),
          "Add a new species")

  .method("choose_cell_in", (List (Simulation::*)(const std::string&))(&Simulation::choose_cell_in),
          "Randomly choose one cell in a genotype")

  .method("choose_cell_in", (List (Simulation::*)(const std::string&,
                                                  const std::vector<RS::AxisPosition>&,
                                                  const std::vector<RS::AxisPosition>&))(&Simulation::choose_cell_in),
          "Randomly choose one cell having a specified genotype in a rectangular selection")

  // schedule_genotype_mutation
  .method("schedule_genotype_mutation", &Simulation::schedule_genotype_mutation,
          "Add a timed mutation between two different species")

  // get_species
  .method("get_species", &Simulation::get_species,
          "Get the species added to the simulation")

  // get_samples_forest
  .method("get_samples_forest", &Simulation::get_samples_forest,
          "Get the descendants forest having as leaves the sampled cells")

  // death_activation_level
  .property("death_activation_level", &Simulation::get_death_activation_level,
                                      &Simulation::set_death_activation_level,
            "The number of cells in a species that activates cell death" )

  // duplicate_internal_cells
  .property("duplicate_internal_cells", &Simulation::get_duplicate_internal_cells,
                                        &Simulation::set_duplicate_internal_cells,
            "Enable/disable duplication for internal cells" )

  // get_clock
  .method("get_clock", &Simulation::get_clock, "Get the current simulation time")

  // get_cell
  .method("get_cell", (List (Simulation::*)(const RS::AxisPosition&,
                                            const RS::AxisPosition&) const)(&Simulation::get_cell),
          "Get one cell from the simulated tissue")

  // get_cells
  .method("get_cells", (List (Simulation::*)(const std::vector<RS::AxisPosition>&,
                                             const std::vector<RS::AxisPosition>&,
                                             const std::vector<std::string>&,
                                             const std::vector<std::string>&) const)(&Simulation::get_cells),
          "Get cells from the simulated tissue")
  .method("get_cells", (List (Simulation::*)(const SEXP&, const SEXP&) const)(&Simulation::get_cells),
          "Get cells from the simulated tissue")
  .method("get_cells", (List (Simulation::*)() const)(&Simulation::get_cells),
          "Get cells from the simulated tissue")

  // get_name
  .method("get_name", &Simulation::get_name, "Get the simulation name")

  // get_lineage_graph
  .method("get_lineage_graph", &Simulation::get_lineage_graph,
          "Get the simulation lineage graph")

  // get_tissue_name
  .method("get_tissue_name", &Simulation::get_tissue_name, "Get the simulation tissue name")

  // get_tissue_size
  .method("get_tissue_size", &Simulation::get_tissue_size, "Get the simulation tissue size")

  // get_added_cells
  .method("get_added_cells", &Simulation::get_added_cells,
          "Get the cells manually added to the simulation")

  // get_counts
  .method("get_counts", &Simulation::get_counts, "Get the current number of cells per species")

  // get_count_history
  .method("get_count_history", (List (Simulation::*)() const)&Simulation::get_count_history,
          "Get the number of simulated events per species along the computation")

  // get_firings
  .method("get_firings", &Simulation::get_firings,
          "Get the current number of simulated events per species")

  // get_firing_history
  .method("get_firing_history", (List (Simulation::*)() const)&Simulation::get_firing_history,
          "Get the number of simulated events per species along the computation")

  // get_rates
  .method("get_rates", &Simulation::get_rates,
          "Get the rates of a species")

  // get_samples_info
  .method("get_samples_info", &Simulation::get_samples_info,
          "Get some pieces of information about the collected samples")

  // history_delta
  .property("history_delta", &Simulation::get_history_delta,
                             &Simulation::set_history_delta,
            "The sampling delta for the get_*_history functions" )

  // mutate
  .method("mutate_progeny",  (void (Simulation::*)(const List&, const std::string&))
                                                  (&Simulation::mutate_progeny),
          "Duplicate a cell and mutate one of its children")
  .method("mutate_progeny",  (void (Simulation::*)(const RS::AxisPosition&,
                                                   const RS::AxisPosition&, const std::string&))
                                                  (&Simulation::mutate_progeny),
          "Duplicate a cell and mutate one of its children")

  // run_up_to_time
  .method("run_up_to_time", &Simulation::run_up_to_time,
          "Simulate the system up to the specified simulation time")

  // run_up_to_event
  .method("run_up_to_event", &Simulation::run_up_to_event,
          "Simulate the system up to the specified number of events")

  // run_up_to_size
  .method("run_up_to_size", &Simulation::run_up_to_size,
          "Simulate the system up to the specified number of cells in the species")

  // sample_cells
  .method("sample_cells", &Simulation::sample_cells,
          "Sample a rectangular region of the tissue")

  // update rates
  .method("update_rates", &Simulation::update_rates,
          "Update the rates of a species")

  // update_tissue
  .method("update_tissue", (void (Simulation::*)(const std::string&, const RS::AxisSize&,
                                                 const RS::AxisSize&))(&Simulation::update_tissue),
          "Update tissue name and size")
  .method("update_tissue", (void (Simulation::*)(const RS::AxisSize&,
                                                 const RS::AxisSize&))(&Simulation::update_tissue),
          "Update tissue size")
  .method("search_sample", &Simulation::search_sample, 
          "Search a rectangular sample containing a given number of cells");

  // recover_simulation
  function("recover_simulation", &Simulation::load,
           "Recover a simulation");

  class_<SamplesForest>("SamplesForest")
    // get_nodes
    .method("get_nodes", (List (SamplesForest::*)() const)(&SamplesForest::get_nodes),
            "Get the nodes of the forest")

    // get_coalescent_cells
    .method("get_coalescent_cells",
            (List (SamplesForest::*)(const std::list<Races::Drivers::CellId>&) const)
                (&SamplesForest::get_coalescent_cells),
            "Get the most recent common ancestor of some cells")

    // get_coalescent_cells
    .method("get_coalescent_cells",
            (List (SamplesForest::*)() const)(&SamplesForest::get_coalescent_cells),
            "Get the most recent common ancestor of all the forest trees")

    // get_subforest_for
    .method("get_subforest_for", &SamplesForest::get_subforest_for,
            "Get the sub-forest for some of the original samples")

    // get_samples_info
    .method("get_samples_info", &SamplesForest::get_samples_info,
            "Get some pieces of information about the samples")

    // get_species
    .method("get_species_info", &SamplesForest::get_species_info,
            "Get the recorded species")

    // show
    .method("show", &SamplesForest::show,
            "Describe the SampleForest");
}

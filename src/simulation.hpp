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

#ifndef __RRACES_SIMULATION__
#define __RRACES_SIMULATION__

#include <vector>
#include <set>
#include <string>

#include <Rcpp.h>

#include <simulation.hpp>

#include "tissue_rectangle.hpp"
#include "samples_forest.hpp"


struct PlainChooser
{
  std::shared_ptr<Races::Drivers::Simulation::Simulation> sim_ptr;
  std::string genotype_name;

  PlainChooser(const std::shared_ptr<Races::Drivers::Simulation::Simulation>& sim_ptr,
               const std::string& genotype_name);

  inline const Races::Drivers::Simulation::CellInTissue& operator()()
  {
    return sim_ptr->choose_cell_in(genotype_name,
                                   Races::Drivers::CellEventType::DUPLICATION);
  }
};

struct RectangularChooser : public PlainChooser
{
  Races::Drivers::RectangleSet rectangle;

  RectangularChooser(const std::shared_ptr<Races::Drivers::Simulation::Simulation>& sim_ptr,
                     const std::string& genotype_name,
                     const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
                     const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner);

  inline const Races::Drivers::Simulation::CellInTissue& operator()()
  {
    return sim_ptr->choose_cell_in(genotype_name, rectangle,
                                   Races::Drivers::CellEventType::DUPLICATION);
  }
};

//' @name Simulation
//' @title Simulates the cell evolution on a tissue
//' @description The objects of this class can simulate the evolution
//'   of many cells belonging to different *species* on a tissue. Each
//'   cell can duplicate or die according to the rates that delineate
//'   the cell species.
//'
//'   `Simulation` supports epigenetic evolutions, and it lets users
//'   define species pairs that have the same genotype (even though,
//'   its genomic characterization is unknown) and differ because
//'   of their epigenetic state (i.e., either "+" or "-").
//'
//'   `Simulation` models epigenetic mutations and allows a cell in
//'   one of a genotype species to generate a new cell belonging to
//'   the other species of the same genotype at a specified rate.
//'
//'   `Simulation` also allows users to schedule mutations from one
//'   genotype to a different genotype.
//' @field add_genotype Adds a genotype and its species \itemize{
//' \item \emph{Parameter:} \code{genotype} - The genotype name.
//' \item \emph{Parameter:} \code{epigenetic_rates} - The epigenetic rates of the genotype species (optional).
//' \item \emph{Parameter:} \code{growth_rates} - The duplication rates of the genotype species.
//' \item \emph{Parameter:} \code{death_rates} - The death rates of the genotype species.
//' }
//' @field choose_cell_in Chooses one cell in a genotype \itemize{
//' \item \emph{Parameter:} \code{genotype} - The genotype of the cell to choose.
//' \item \emph{Parameter:} \code{lower_corner} - The lower left corner of a rectangular selection (optional).
//' \item \emph{Parameter:} \code{upper_corner} - The upper right corner of a rectangular selection (optional).
//' \item \emph{Returns:} A list reporting "cell_id", "genotype", "epistate", "position_x",
//'    and "position_y" of the choosen cell.
//' }
//' @field death_activation_level The number of cells that activates cell death in a species.
//' @field duplicate_internal_cells Enable/disable duplication for internal cells.
//' @field get_added_cells Gets the cells manually added to the simulation \itemize{
//' \item \emph{Returns:} A data frame reporting "genotype", "epistate", "position_x",
//'         "position_y", and "time" for each cells manually added to
//'         the simulation.
//' }
//' @field search_sample Seach a rectangular sample having a minimum number of cells\itemize{
//' \item \emph{Parameter:} \code{genotype_name} - The genotype of the searched cells.
//' \item \emph{Parameter:} \code{num_of_cells} - The number of cells in the searched sample.
//' \item \emph{Parameter:} \code{width} - The width of the searched sample.
//' \item \emph{Parameter:} \code{height} - The height of the searched sample.
//' \item \emph{Returns:} If a rectangular sample satisfying the provided constraints can 
//'               be found, the corresponding rectangle.
//' }
//' @field get_cell Gets one the tissue cells \itemize{
//' \item \emph{Parameter:} \code{x} - The position of the aimed cell on the x axis.
//' \item \emph{Parameter:} \code{y} - The position of the aimed cell on the y axis.
//' \item \emph{Returns:} A data frame reporting "cell_id", "genotype", "epistate", "position_x",
//'    and "position_y" of the aimed cell.
//' }
//' @field get_cells Gets the tissue cells \itemize{
//' \item \emph{Parameter:} \code{lower_corner} - The lower-left corner of the selection frame (optional).
//' \item \emph{Parameter:} \code{upper_corner} - The upper-right corner of the selection frame (optional).
//' \item \emph{Parameter:} \code{genotype_filter} - The vector of the to-be-selected genotype names (optional).
//' \item \emph{Parameter:} \code{epigenetic_filter} - The vector of the to-be-selected epigenetic states (optional).
//' \item \emph{Returns:} A data frame reporting "cell_id", "genotype", "epistate", "position_x",
//'    and "position_y" for each cells satisfying the provided filters and laying
//'    in the input frame.
//' }
//' @field get_clock Gets the simulated time \itemize{
//' \item \emph{Returns:} The time simulated by the simulation.
//' }
//' @field get_count_history Gets the history of the number of cells per species \itemize{
//' \item \emph{Returns:} A data frame reporting "genotype", "epistate", "counts",
//'     and "time" for each species and for each sampled time.
//' }
//' @field get_counts Counts the number of cells \itemize{
//' \item \emph{Returns:} A data frame reporting "genotype", "epistate", "counts" for each
//'      species in the simulation.
//' }
//' @field get_firing_history Gets the history of the number of fired events \itemize{
//' \item \emph{Returns:} A data frame reporting "event", "genotype", "epistate", "fired",
//'      and "time" for each event type, for each species, and for each sampled time.
//' }
//' @field get_firings Gets the number of fired events \itemize{
//' \item \emph{Returns:} A data frame reporting "event", "genotype", "epistate", and "fired"
//'     for each event type and for each species.
//' }
//' @field get_name Gets the simulation name \itemize{
//' \item \emph{Returns:} The simulation name, which corresponds to the name of the directory
//'         in which the simulation is saving its progresses.
//' }
//' @field get_lineage_graph Gets the simulation lineage graph\itemize{
//' \item \emph{Returns:} A data frame reporting "ancestor", "progeny", and "first_occurrence"
//'         of each species-to-species transition.
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
//'         during the simulation, the columns "name", "time", "ymin",
//'         "xmin", "ymax", "xmax", and  "tumoral cells". "ymin",
//'         "xmin", "ymax", "xmax" report the boundaries of the sampled
//'         rectangular region, while "tumoral cells" is the number of
//'         tumoral cells in the sample.
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
//' \item \emph{Parameter:} \code{mutated_genotype} - The genotype of the mutated cell.
//' }
//' or
//' \itemize{
//' \item \emph{Parameter:} \code{x} - The position of the cell whose progeny will mutate on the x axis.
//' \item \emph{Parameter:} \code{y} - The position of the cell whose progeny will mutate on the y axis.
//' \item \emph{Parameter:} \code{mutated_genotype} - The genotype of the mutated cell.
//' }
//' @field place_cell Place one cell in the tissue \itemize{
//' \item \emph{Parameter:} \code{species} - The name of the new cell species.
//' \item \emph{Parameter:} \code{x} - The position on the x axis of the cell.
//' \item \emph{Parameter:} \code{y} - The position on the y axis of the cell.
//' }
//' @field schedule_genotype_mutation Schedules a genotype mutation \itemize{
//' \item \emph{Parameter:} \code{src} - The name of the genotype from which the mutation occurs.
//' \item \emph{Parameter:} \code{dest} - The name of the genotype to which the mutation leads.
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
class Simulation
{
  std::shared_ptr<Races::Drivers::Simulation::Simulation> sim_ptr;  //!< The pointer to a RACES simulation object
  std::string name;      //!< The simulation name
  bool save_snapshots;   //!< A flag to preserve binary dump after object destruction

  void init(const SEXP& sexp);

  static bool has_names(const Rcpp::List& list, std::vector<std::string> aimed_names);

  static bool has_names_in(const Rcpp::List& list, std::set<std::string> aimed_names);

  static std::vector<Races::Drivers::Simulation::Direction> get_possible_directions();

  Rcpp::List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                 const std::set<Races::Drivers::SpeciesId> &species_filter,
                 const std::set<std::string> &epigenetic_filter) const;

  Rcpp::List wrap_a_cell(const Races::Drivers::Simulation::CellInTissue& cell) const;

public:
  Simulation();

  Simulation(const SEXP& sexp);

  Simulation(const SEXP& first_param, const SEXP& second_param);

  Simulation(const std::string& simulation_name, const int& seed, const bool& save_snapshots);

  ~Simulation();

  inline void update_tissue(const std::string& name,
                            const Races::Drivers::Simulation::AxisSize& width,
                            const Races::Drivers::Simulation::AxisSize& height)
  {
    sim_ptr->set_tissue(name, {width, height});
  }

  inline void update_tissue(const Races::Drivers::Simulation::AxisSize& width,
                            const Races::Drivers::Simulation::AxisSize& height)
  {
    sim_ptr->set_tissue("A tissue", {width, height});
  }

  void add_genotype(const std::string& genotype, const Rcpp::List& epigenetic_rates,
                    const Rcpp::List& growth_rates, const Rcpp::List& death_rates);

  void add_genotype(const std::string& genotype, const double& growth_rate, const double& death_rate);

  Rcpp::List add_genotype(const Rcpp::List& prova);

  inline Races::Time get_clock() const
  {
    return sim_ptr->get_time();
  }

  void place_cell(const std::string& species_name,
                  const Races::Drivers::Simulation::AxisPosition& x,
                  const Races::Drivers::Simulation::AxisPosition& y);

  size_t count_history_sample_in(const Races::Time& minimum_time,
                                 const Races::Time& maximum_time) const;

  Rcpp::List get_added_cells() const;

  Rcpp::List get_counts() const;

  inline Rcpp::List get_count_history() const
  {
    return get_count_history(0);
  }

  Rcpp::List get_count_history(const Races::Time& minimum_time) const;

  Rcpp::List get_count_history(const Races::Time& minimum_time,
                         const Races::Time& maximum_time) const;

  Rcpp::List get_cells() const;

  Rcpp::List get_cell(const Races::Drivers::Simulation::AxisPosition& x,
                const Races::Drivers::Simulation::AxisPosition& y) const;

  Rcpp::List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner) const;

  Rcpp::List get_cells(const SEXP& first_param, const SEXP& second_param) const;

  Rcpp::List get_cells(const std::vector<std::string>& species_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  Rcpp::List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                 const std::vector<std::string>& genotype_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  Rcpp::List get_lineage_graph() const;

  Rcpp::List get_samples_info() const;

  inline void schedule_genotype_mutation(const std::string& src, const std::string& dst,
                                         const Races::Time& time)
  {
    sim_ptr->schedule_genotype_mutation(src, dst, time);
  }

  void run_up_to_time(const Races::Time& time);

  void run_up_to_size(const std::string& species_name, const size_t& num_of_cells);

  void run_up_to_event(const std::string& event, const std::string& species_name,
                       const size_t& num_of_events);

  void sample_cells(const std::string& sample_name,
                    const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
                    const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner) const;

  Rcpp::List get_firings() const;

  inline Rcpp::List get_firing_history() const
  {
    return get_firing_history(0);
  }

  Rcpp::List get_firing_history(const Races::Time& minimum_time) const;

  Rcpp::List get_firing_history(const Races::Time& minimum_time,
                          const Races::Time& maximum_time) const;

  Rcpp::List get_species() const;

  inline std::string get_name() const
  {
    return name;
  }

  inline const std::string& get_tissue_name() const
  {
    return sim_ptr->tissue().get_name();
  }

  Rcpp::IntegerVector get_tissue_size() const;

  Rcpp::List get_rates(const std::string& species_name) const;

  void update_rates(const std::string& species_name, const Rcpp::List& list);

  template<typename CHOOSER, std::enable_if_t<std::is_base_of_v<PlainChooser, CHOOSER>, bool> = true>
  Rcpp::List choose_border_cell_in(CHOOSER& chooser)
  {
    namespace RS = Races::Drivers::Simulation;

    const auto directions = get_possible_directions();

    const RS::Tissue& tissue = sim_ptr->tissue();

    size_t i{0};
    while (++i<1000) {
      const auto& cell = chooser();

      for (const auto& dir: directions) {
        RS::PositionInTissue pos = cell;

        do {
          pos = pos + RS::PositionDelta(dir);
        } while (tissue.is_valid(pos) && tissue(pos).is_wild_type());

        if (!tissue.is_valid(pos)) {
          return wrap_a_cell(cell);
        }
      }
    }

    throw std::domain_error("Missed to find a border cell");
  }

  Rcpp::List choose_border_cell_in(const std::string& genotype_name);

  Rcpp::List choose_border_cell_in(const std::string& genotype_name,
                             const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
                             const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner);

  Rcpp::List choose_cell_in(const std::string& genotype_name);

  Rcpp::List choose_cell_in(const std::string& genotype_name,
                      const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
                      const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner);

  void mutate_progeny(const Races::Drivers::Simulation::AxisPosition& x,
                      const Races::Drivers::Simulation::AxisPosition& y,
                      const std::string& mutated_genotype);

  void mutate_progeny(const Rcpp::List& cell_position, const std::string& mutated_genotype);

  inline size_t get_death_activation_level() const
  {
    return sim_ptr->death_activation_level;
  }

  inline void set_death_activation_level(const size_t death_activation_level)
  {
    sim_ptr->death_activation_level = death_activation_level;
  }

  inline bool get_duplicate_internal_cells() const
  {
    return sim_ptr->duplicate_internal_cells;
  }

  inline void set_duplicate_internal_cells(const bool duplicate_internal_cells)
  {
    sim_ptr->duplicate_internal_cells = duplicate_internal_cells;
  }

  inline Races::Time get_history_delta() const
  {
    return sim_ptr->get_statistics().get_history_delta();
  }

  inline void set_history_delta(const Races::Time history_time_delta)
  {
    sim_ptr->get_statistics().set_history_delta(history_time_delta);
  }

  static Simulation load(const std::string& directory_name);

  SamplesForest get_samples_forest() const;

  TissueRectangle search_sample(const std::string& genotype_name, const size_t& num_of_cells,
                                const uint16_t& width, const uint16_t& height);
};

//' @name Simulation$new
//' @title Constructs a new Simulation
//' @param simulation_name The name of the simulation (optional).
//' @param seed The seed for the pseudo-random generator (optional).
//' @param save_snapshots A flag to save simulation snapshots on disk (optional, 
//'                default `FALSE`).
//' @examples
//' # create a Simulation object storing binary dump in a temporary directory. 
//' # The data are deleted from the disk as soon as the object is destroyed.
//' sim <- new(Simulation, "test")
//'
//' # add a new species, place a cell in the tissue, and let the simulation evolve.
//' sim$add_genotype(genotype = "A", growth_rate = 0.3, death_rate = 0.02)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_time(30)
//'
//' # no directory "test" has been created
//' "test" %in% list.files(".")
//'
//' # (let us delete the directory "test" manually)
//' unlink("test", recursive = TRUE)
//' 
//' # By using the optional parameter `save_snapshots`, we force the 
//' # simulation to save its progresses in a local directory whose name 
//' # is the name of the simulation, i.e., "test". This data will be 
//' # preserved when the simulation object will be destroyed.
//' sim <- new(Simulation, "test", save_snapshots=TRUE)
//'
//' # as done above, we add a new species, place a cell in the tissue, and let the 
//' # simulation evolve.
//' sim$add_genotype(genotype = "A", growth_rate = 0.3, death_rate = 0.02)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_time(30)
//'
//' # the directory "test" exists and contains a binary dump of 
//' # sthe simulation.
//' "test" %in% list.files(".")
//'
//' # let us manually delete the "test" directory
//' unlink("test", recursive=TRUE)
//'
//' # we can also provide a random seed to the simulation...
//' sim <- new(Simulation, "test", 13)
//'
//' # ...or creating a simulation without providing any name. By default, the 
//' # simulation name will have the following format `races_<date>_<hour>`.
//' sim <- new(Simulation, 13)

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

//' @name Simulation$add_genotype
//' @title Adds a genotype and its species
//' @description This method adds a genotype and its species to the
//'      simulation. If the optional parameter `epigenetic_rate` is
//'      provided, then two new species having the same genotype and
//'      opposite epigenetic states are created. When, instead, the
//'      optional parameter `epigenetic_rate` is missing, this
//'      method creates only one species with no epigenetic states.
//' @param genotype The genotype name.
//' @param epigenetic_rates The epigenetic rates of the genotype species (optional).
//' @param growth_rates The duplication rates of the genotype species.
//' @param death_rates The death rates of the genotype species.
//' @examples
//' sim <- new(Simulation)
//'
//' # create the two species "A+" and "A-". They both have genotype "A".
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # create the species "C" its genotype is "C".
//' sim$add_genotype(genotype = "C", growth_rate = 0.2, death_rate = 0.1)


//' @name Simulation$get_species
//' @title Gets the species
//' @return A data frame reporting "genotype", "epistate", "growth_rate",
//'    "death_rate", and "switch_rate" for each registered species.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_genotype("B", growth_rate = 0.15, death_rate = 0.05)
//'
//' # get the added species and their rates. In this case, "A"
//' # and "B"
//' sim$get_species()


//' @name Simulation$place_cell
//' @title Place one cell in the tissue
//' @param species The name of the new cell species.
//' @param x The position on the x axis of the cell.
//' @param y The position on the y axis of the cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # add into the tissue a cell of species "A+" in position (500,500)
//' sim$place_cell("A+", 500, 500)


//' @name Simulation$get_cell
//' @title Gets one of the tissue cells
//' @description This method collects some data of the aimed cell without altering
//'      the tissue.
//' @param x The position of the aimed cell on the x axis.
//' @param y The position of the aimed cell on the y axis.
//' @return A data frame reporting "cell_id", "genotype", "epistate", "position_x",
//'    and "position_y" of the aimed cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.02, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.3, "-" = 0.1),
//'                  death_rates = c("+" = 0.02, "-" = 0.01))
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(40)
//'
//' # collect all the cells in the tissue
//' sim$get_cell(501, 502)

  
//' @name Simulation$get_cells
//' @title Gets the tissue cells
//' @description This method collects some data about the cells in the tissue
//'      without altering the tissue itself. The pairs of optional parameters
//'      `lower_corner` and `upper_corner` define a frame of the tissue in
//'      which the data are sampled. The optional parameters `genotype_filter`
//'      and `epigenetic_filter` filter the collected cell data according to
//'      the cell genotype and epigenetic state.
//' @param lower_corner The lower-left corner of the selection frame (optional).
//' @param upper_corner The upper-right corner of the selection frame (optional).
//' @param genotype_filter The vector of the to-be-selected genotype names (optional).
//' @param epigenetic_filter The vector of the to-be-selected epigenetic states (optional).
//' @return A data frame reporting "cell_id", "genotype", "epistate", "position_x",
//'    and "position_y" for each cells satisfying the provided filters and laying
//'    in the input frame.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.3, "-" = 0.1),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(30)
//'
//' # collect all the cells in the tissue
//' sim$get_cells()
//'
//' # get the cells in the frame [495,505]x[490,500]
//' sim$get_cells(lower_corner=c(495,490), upper_corner=c(505,500))
//'
//' # cells can be filtered by genotype name...
//' sim$get_cells(genotype_filter=c("A"),epigenetic_filter=c("+","-"))
//'
//' # ...or by epigenetic state
//' sim$get_cells(genotype_filter=c("A","B"),epigenetic_filter=c("-"))
//'
//' # cells can be filtered by frame, genotype, and epigenetic states
//' sim$get_cells(lower_corner=c(495,495), upper_corner=c(505,505),
//'               genotype_filter=c("A"),epigenetic_filter=c("+","-"))


//' @name Simulation$get_counts
//' @title Counts the number of cells
//' @return A data frame reporting "genotype", "epistate", "counts" for each
//'      species in the simulation.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_genotype("B", growth_rate = 0.15, death_rate = 0.05)
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_counts()


//' @name Simulation$get_added_cells
//' @title Gets the cells manually added to the simulation
//' @return A data frame reporting "genotype", "epistate", "position_x",
//'         "position_y", and "time" for each cells manually added to
//'         the simulation.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.3, "-" = 0.1),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 30)
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(50)
//'
//' # counts the number of cells per species
//' sim$get_added_cells()


//' @name Simulation$schedule_genotype_mutation
//' @title Schedules a genotype mutation
//' @description This method schedules a genotype mutation that can occur
//'      from any of the species of the source genotype to the species of
//'      the destination genotype with a consistent epigenetic state.
//'      For the sake of example, if the mutation from "A" to "B" is
//'      scheduled, then we have three possible situations:
//'      1. The genotype "A" consists of the only species "A". Then,
//'         during one duplication of a cell of "A", one cell of "B"
//'         will arise.
//'      2. The genotype "A" consists of the species "A+" and "A-" and
//'         during one duplication of a cell of "A+", one cell of "B+"
//'         will arise.
//'      3. The genotype "A" consists of the species "A+" and "A-" and
//'         during one duplication of a cell of "A-", one cell of "B-"
//'         will arise.
//'      No other scenario can occur.
//' @param src The name of the genotype from which the mutation occurs.
//' @param dest The name of the genotype to which the mutation leads.
//' @param time The simulated time at which the mutation will occurs.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.3, "-" = 0.1),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # schedule an evolution from genotype "A" to genotype "B" at time 50
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)


//' @name Simulation$get_lineage_graph
//' @title Gets the simulation lineage graph
//' @description At the beginning of the computation only the species of the added
//'         cells are present in the tissue. As the simulation proceeds new species
//'         arise as a consequence of either genotype mutations or epigenetic
//'         switches. The *lineage graph* stores these species evolutions and it
//'         reports the first occurrence time of any species-to-species transition.
//'
//'         This method returns the lineage graph of the simulation.
//' @return A data frame reporting "ancestor", "progeny", and "first_occurrence" of
//'         each species-to-species transition.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.3, "-" = 0.1),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 20)
//' sim$run_up_to_time(50)
//'
//' sim$get_lineage_graph()


//' @name Simulation$run_up_to_time
//' @title Simulates cell evolution
//' @param time The final simulation time.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$place_cell("A", 500, 500)
//'
//' # simulate the tissue up to simulate timed 100
//' sim$run_up_to_time(40)


//' @name Simulation$run_up_to_size
//' @title Simulates cell evolution
//' @description This method simulates cell evolution until the number of cells in
//'       a species reaches a specified threshold.
//' @param species The species whose number of cells is considered.
//' @param num_of_cells The threshold for the cell number.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//'
//' # simulate the tissue until the species "A+" account for 100
//' # contemporary cells
//' sim$run_up_to_size(species = "A+", num_of_cells = 100)


//' @name Simulation$run_up_to_event
//' @title Simulates cell evolution
//' @description This method simulates cell evolution until the number of events that
//'         have occurred to cells of a species reaches a specified threshold.
//' @param event The considered event, i.e., "growth", "death", or "switch".
//' @param species The species whose event number is considered.
//' @param num_of_events The threshold for the event number.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//'
//' # simulate the cell evolution until the number of epigenetic events from
//' # the species "A+" is less than 100.
//' sim$run_up_to_event(event = "switch", species = "A+", num_of_events = 100)


//' @name Simulation$get_clock
//' @title Gets the simulated time
//' @return The time simulated by the simulation.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the simulated time
//' sim$get_clock()


//' @name Simulation$get_firings
//' @title Gets the number of fired events
//' @return A data frame reporting "event", "genotype", "epistate", and "fired"
//'     for each event type, genotype, and epigenetic states.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the number of event fired per event and species
//' sim$get_firings()


//' @name Simulation$get_firing_history
//' @title Gets the history of the number of fired events
//' @description This method returns a data frame reporting the number of
//'           events fired up to each sampled simulation time.
//' @return A data frame reporting "event", "genotype", "epistate", "fired",
//'     and "time" for each event type, for each species, and for each
//'     sampled time.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$history_delta <- 20
//' sim$run_up_to_time(70)
//'
//' # get the number of event fired per event and species
//' sim$get_firing_history()


//' @name Simulation$get_count_history
//' @title Gets the history of the number of cells per species
//' @description This method returns a data frame reporting the number of
//'           species cells in each sampled simulation time.
//' @return A data frame reporting "genotype", "epistate", "counts",
//'     and "time" for each species, and for each sampled time.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_genotype("B", growth_rate = 0.15, death_rate = 0.05)
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)
//' sim$place_cell("A", 500, 500)
//' sim$history_delta <- 20
//' sim$run_up_to_time(70)
//'
//' # get the history of species counts
//' sim$get_count_history()


//' @name Simulation$get_name
//' @title Gets the simulation name
//' @return The simulation name, which corresponds to the name of the directory
//'         in which the simulation is saving its progresses.
//' @examples
//' sim <- new(Simulation)
//'
//' # Expecting "test"
//' sim$get_name()

//' @name Simulation$get_tissue_name
//' @title Gets the tissue name
//' @return The name of the simulated tissue.
//' @examples
//' sim <- new(Simulation)
//' sim$update_tissue("Liver", 1200, 900)
//'
//' # get the tissue name, i.e., expecting "Liver"
//' sim$get_tissue_name()


//' @name Simulation$get_tissue_size
//' @title Gets the size of the simulated tissue
//' @return The vector `c(x_size, y_size)` of the simulated tissue.
//' @examples
//' sim <- new(Simulation)
//' sim$update_tissue("Liver", 1200, 900)
//'
//' # get the tissue size, i.e., expecting c(1200,900)
//' sim$get_tissue_size()


//' @name Simulation$get_rates
//' @title Get the rates of a species
//' @param species The species whose rates are aimed.
//' @return The list of the species rates.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Get the rates of "A-". In this case c("growth"=0.08, "death"=0.01, "switch"=0.02) is expected
//' sim$get_rates("A-")

//' @name Simulation$update_rates
//' @title Update the rates of a species
//' @param species The species whose rates must be updated.
//' @param rates The list of rates to be updated.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Set the death and epigenetic switch rates of "A-" to 0
//' sim$update_rates("A-", c(switch=0, death=0))


//' @name Simulation$choose_cell_in
//' @title Chooses one cell in a genotype
//' @description This method chooses one of the cells whose genotype
//'         is `genotype`. Optionally, the lower and upper corners
//'         of a tissue rectangular selection can be provided
//'         to obtain one cell in the rectangle.
//' @param genotype The genotype of the cell to choose.
//' @param lower_corner The lower corner of the rectangular selection (optional).
//' @param upper_corner The upper corner of the rectangular selection (optional).
//' @return A list reporting "cell_id", "genotype", "epistate", "position_x",
//'    and "position_y" of the choosen cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.1, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.15, "-" = 0.3),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$death_activation_level <- 100
//' sim$schedule_genotype_mutation("A","B",20)
//' sim$run_up_to_size(species = "B-", num_of_cells = 50)
//'
//' # Randomly choose one cell in "B" in the tissue
//' sim$choose_cell_in(genotype = "B")


//' @name Simulation$mutate_progeny
//' @title Generate a mutated progeny
//' @description This method simulates both the duplication of the cell in the
//'       specified position and the birth of one cells of a given
//'       genotype that preserves the epigenetic status of the original cell.
//'       The mutated cell will be located in the position of its parent.
//' @param cell_position The position of the cell whose offspring will mutate.
//' @param mutated_genotype The genotype of the mutated cell.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.01, "-" = 0.01))
//' sim$place_cell("A+", 500, 500)
//' sim$run_up_to_time(30)
//'
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.1, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.15, "-" = 0.3),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # duplicate the cell in position (503, 492). One of
//' # its direct descendents will have genotype "B"
//' # sim$mutate_progeny(503, 492, "B")
//'
//' # the output of `choose_cell_in` and `get_cell` can also be used
//' # as input for `mutate_progeny`
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")


//' @name Simulation$sample_cells
//' @title Sample a tissue rectangle region.
//' @description This method removes a rectangular region from the simulated
//'       tissue and stores its cells in a sample that can subsequently
//'       retrieved to build a samples forest.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  growth_rate = 0.2,
//'                  death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))


//' @name Simulation$get_samples_info
//' @title Retrieve information about the samples
//' @description This method retrieves information about
//'           the samples collected along the simulation.
//'           It returns a data frame reporting, for each
//'           sample, the name, the sampling time, the
//'           position, and the number of tumoural cells.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  growth_rate = 0.2,
//'                  death_rate = 0.01)
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


//' @name Simulation$death_activation_level
//' @title The number of cells that activates cell death in a species.
//' @description This value is the minimum number of cells that
//'       enables cell death in a species. The cell of a species $S$ can die
//'       if and only if that $S$ has reached the death activation level at
//'       least once during the simulation.
//' @examples
//' sim <- new(Simulation)
//'
//' # get the simulation death activation level
//' sim$death_activation_level
//'
//' # set the death activation level to 50
//' sim$death_activation_level <- 50


//' @name Simulation$duplicate_internal_cells
//' @title Enable/disable duplication for internal cells.
//' @description This Boolean flag enable/disable duplication of internal
//'            cells. When it is set to `FALSE`, the border-growth model
//'            is used. Otherwise, the homogeneous-growth model is applied.
//'            It is set to `FALSE` by default.
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

//' @name Simulation$history_delta
//' @title The delta time between time series samples
//' @description This value is the maximum time between two successive
//'          time series data samples.
//' @examples
//' sim <- new(Simulation)
//'
//' # get the delta time between two time series samples (0 by default)
//' sim$history_delta
//'
//' # set the delta time between two time series samples
//' sim$death_activation_level <- 20


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
//' sim$add_genotype("A",
//'                  epigenetic_rates=c("+-" = 0.01, "-+"=0.01),
//'                  growth_rates = c("+"=0.1, "-"=0.01),
//'                  death_rates = c("+"=0.05, "-"=0.005))
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
//' unlink("recover_simulation_test", recursive=TRUE)


//' @name Simulation$get_samples_forest
//' @title Get the samples forest
//' @return The samples forest having as leaves the sampled cells
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  growth_rate = 0.2,
//'                  death_rate = 0.01)
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


//' @name Simulation$search_sample
//' @title Search a rectangular sample containing a minimum number of cells
//' @description This method searches a rectangular tissue sample containing 
//'        the provided number of cells. The sizes of the sample are also
//'        provided a parameter of the method. 
//'        The complexity of this method is O(|tissue rows|*|tissue cols|).
//' @param genotype_name The genotype of the searched cells.
//' @param num_of_cells The number of cells in the searched sample.
//' @param width The width of the searched sample.
//' @param height The height of the searched sample.
//' @return If a rectangular sample satisfying the provided constraints can 
//'               be found, the corresponding rectangle.
//' @examples
//' sim <- new(Simulation)
//' sim$death_activation_level <- 50
//' sim$add_genotype(genotype = "A", growth_rate = 0.2, death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//' sim$run_up_to_size(species = "A", num_of_cells = 500)
//'
//' sim$add_genotype(genotype = "B", growth_rate = 0.3, death_rate = 0.01)
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
//' sim$run_up_to_size(species = "B", num_of_cells = 1000)
//'
//' # find a 10x10 sample containing 80 "B" cells
//' sim$search_sample("B",80,50,50)

#endif // __RRACES_SIMULATION__

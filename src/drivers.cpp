/* 
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023 Alberto Casagrande.
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
#include <set>

#include <Rcpp.h>

#include "simulation.hpp"
#include "ending_conditions.hpp"

using namespace Rcpp;

template<typename SIMULATION_TEST>
struct RTest : public SIMULATION_TEST
{
  size_t counter;

  template<typename ...Args>
  explicit RTest(Args...args):
      SIMULATION_TEST(args...), counter{0}
  {} 

  bool operator()(const Races::Drivers::Simulation::Simulation& simulation)
  {
    if (++counter >= 10000) {
      counter = 0;
      try {
        Rcpp::checkUserInterrupt();
      } catch (...) {
        return true;
      }
    }

    using namespace Races::Drivers::Simulation;

    return SIMULATION_TEST::operator()(simulation);
  }
};

const std::map<std::string, Races::Drivers::CellEventType> event_names{
  {"death",  Races::Drivers::CellEventType::DEATH},
  {"growth", Races::Drivers::CellEventType::DUPLICATION},
  {"switch", Races::Drivers::CellEventType::EPIGENETIC_SWITCH},
};

size_t count_events(const Races::Drivers::Simulation::SpeciesStatistics& statistics,
                    const Races::Drivers::CellEventType& event)
{
  switch(event)
  {
    case Races::Drivers::CellEventType::DEATH:
      return statistics.killed_cells;
    case Races::Drivers::CellEventType::DUPLICATION:
      return statistics.num_duplications;
    case Races::Drivers::CellEventType::EPIGENETIC_SWITCH:
      return statistics.num_of_epigenetic_events();
    default:
      throw std::domain_error("get_counts: unsupported event");
  }
}

void handle_unknown_event(const std::string& event)
{
  std::ostringstream oss;

  oss << "Event \"" << event << "\" is not supported. " << std::endl
      << "Supported events are ";

  size_t i{0};
  for (const auto& [name, type] : event_names) {
    if (i>0) {
      if (event_names.size()!=2) {
        oss << ",";
      }
    }

    if ((++i)+1==event_names.size()) {
      oss << " and ";
    }

    oss << "\"" << name << "\"";
  }

  oss << ".";

  throw std::domain_error(oss.str());
}

std::set<Races::Drivers::SpeciesId> 
get_species_ids_from_genotype_name(const Races::Drivers::Simulation::Tissue& tissue,
                                   const std::set<std::string>& genotype_name)
{
  std::set<Races::Drivers::SpeciesId > species_ids;

  for (const auto& species: tissue) {
    if (genotype_name.count(species.get_genotype_name())>0) {
      species_ids.insert(species.get_id());
    }
  }

  return species_ids;
}

Races::Drivers::Simulation::PositionInTissue 
get_position_in_tissue(const std::vector<Races::Drivers::Simulation::AxisPosition>& position)
{
  if (position.size()==2) {
    return {position[0], position[1]};
  }

  throw std::domain_error("rRACES supports only 2 dimensional space so far");
}

Races::Drivers::RectangleSet 
get_rectangle(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner,
              const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner)
{
  auto l_position = get_position_in_tissue(lower_corner);
  auto u_position = get_position_in_tissue(upper_corner);

  return {l_position, u_position};
}

size_t count_driver_mutated_cells(const Races::Drivers::Simulation::Tissue& tissue,
                                  const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                                  const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                                  const std::set<Races::Drivers::SpeciesId>& species_filter,
                                  const std::set<std::string>& epigenetic_filter)
{
  using namespace Races::Drivers;
  using namespace Races::Drivers::Simulation;

  if (lower_corner.size() != upper_corner.size()) {
    throw std::runtime_error("lower_corner and upper_corner must have the same size");
  }
  
  auto lower_it = lower_corner.begin();
  auto upper_it = upper_corner.begin();
  for (;lower_it!=lower_corner.end();++lower_it,++upper_it) 
  {
    if (*lower_it>*upper_it) {
      return 0;
    }
  }

  size_t total{0};
  for (auto x=lower_corner[0]; x<=upper_corner[0]; ++x) {
    for (auto y=lower_corner[1]; y<=upper_corner[1]; ++y) {
      auto cell_proxy = tissue({x,y});
      if (!cell_proxy.is_wild_type()) {
        const CellInTissue& cell = cell_proxy;

        if (species_filter.count(cell.get_species_id())>0) {

          const auto& species = tissue.get_species(cell.get_species_id());
          const auto& signature = species.get_methylation_signature();
          auto sign_string = GenotypeProperties::signature_to_string(signature);

          if (epigenetic_filter.count(sign_string)>0) {
            ++total;
          }
        }
      }
    }
  }

  return total;
}

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
//' @field add_cell Implant one cell in the tissue \itemize{
//' \item \emph{Parameter:} \code{species} - The name of the new cell species.
//' \item \emph{Parameter:} \code{x} - The position on the x axis of the cell.
//' \item \emph{Parameter:} \code{y} - The position on the y axis of the cell.
//' }
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
//' @field get_counts Counts the number of cells \itemize{
//' \item \emph{Returns:} A data frame reporting "genotype", "epistate", "counts" for each
//'      species in the simulation.
//' }
//' @field get_directory Gets the simulation directory \itemize{
//' \item \emph{Returns:} The directory in which the simulation is saving its progress.
//' }
//' @field get_firings Gets the number of fired events \itemize{
//' \item \emph{Returns:} A data frame reporting "event", "genotype", "epistate", and "fired"
//'     for each event type, genotype, and epigenetic states.
//' }
//' @field get_rates Gets the rates of a species\itemize{
//' \item \emph{Parameter:} \code{species} - The species whose rates are aimed.
//' \item \emph{Returns:} The list of the species names.
//' }
//' @field get_species_names Gets the species name \itemize{
//' \item \emph{Returns:} The vector of the species names.
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
class Simulation : private Races::Drivers::Simulation::Simulation
{
  static bool has_names(const List& list, std::vector<std::string> aimed_names)
  {
    if (aimed_names.size() != static_cast<size_t>(list.size())) {
      return false;
    }

    for (const std::string& name: aimed_names) {
      if (!list.containsElementNamed(name.c_str())) {
        return false;
      }
    }

    return true;
  }

  static bool has_names_in(const List& list, std::set<std::string> aimed_names)
  {
    if (aimed_names.size() < static_cast<size_t>(list.size())) {
      return false;
    }

    CharacterVector names = wrap(list.names());

    for (size_t i=0; i<static_cast<size_t>(names.size()); ++i) {
      if (aimed_names.count(as<std::string>(names[i]))==0) {
        return false;
      }
    }

    return true;
  }

  List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                 const std::set<Races::Drivers::SpeciesId> &species_filter,
                 const std::set<std::string> &epigenetic_filter) const
  {
    namespace RS = Races::Drivers::Simulation;

    using namespace Races::Drivers;

    if (lower_corner.size() != 2) {
      throw std::domain_error("The lower corner must be a vector having size 2");
    }

    if (upper_corner.size() != 2) {
      throw std::domain_error("The upper corner must be a vector having size 2");
    }

    size_t num_of_rows = count_driver_mutated_cells(tissue(), lower_corner, upper_corner,
                                                    species_filter, epigenetic_filter);

    IntegerVector ids(num_of_rows);
    CharacterVector genotype_names(num_of_rows);
    CharacterVector epi_stati(num_of_rows);
    IntegerVector x_pos(num_of_rows);
    IntegerVector y_pos(num_of_rows);

    size_t i{0};
    for (auto x=lower_corner[0]; x<=upper_corner[0]; ++x) {
      for (auto y=lower_corner[1]; y<=upper_corner[1]; ++y) {
        auto cell_proxy = tissue()({x,y});
        if(!cell_proxy.is_wild_type()) {

          const RS::CellInTissue& cell = cell_proxy;

          const auto& species = tissue().get_species(cell.get_species_id());
          const auto& signature = species.get_methylation_signature();
          auto sign_string = GenotypeProperties::signature_to_string(signature);
          
          if (species_filter.count(cell.get_species_id())>0
               && epigenetic_filter.count(sign_string)>0) { 

            ids[i] = cell.get_id();
            genotype_names[i] = species.get_genotype_name();

            epi_stati[i] = GenotypeProperties::signature_to_string(signature);

            x_pos[i] = x;
            y_pos[i] = y;

            ++i;
          }
        }
      }
    }

    return DataFrame::create(_["cell_id"]=ids, _["genotype"]=genotype_names,
                             _["epistate"]=epi_stati, _["position_x"]=x_pos,
                             _["position_y"]=y_pos);
  }

  List wrap_a_cell(const Races::Drivers::Simulation::CellInTissue& cell) const
  {
    using namespace Races::Drivers;

    const auto& species = tissue().get_species(cell.get_species_id());

    const auto& genotype_name = find_genotype_name(species.get_genotype_id());

    auto epistate = GenotypeProperties::signature_to_string(species.get_methylation_signature());

    return List::create(_["cell_id"]=cell.get_id(), _["genotype"]=genotype_name,
                        _["epistate"]=epistate, _["position_x"]=cell.x,
                        _["position_y"]=cell.y);
  }

public:
  Simulation();

  Simulation(const int& seed);

  Simulation(const std::string& output_dir);

  Simulation(const std::string& output_dir, const int& seed);

  void update_tissue(const std::string& name, const uint16_t& width, const uint16_t& height);

  void update_tissue(const uint16_t& width, const uint16_t& height);

  void add_genotype(const std::string& genotype, const List& epigenetic_rates, 
                    const List& growth_rates, const List& death_rates);

  void add_genotype(const std::string& genotype, const double& growth_rate, const double& death_rate);

  inline Races::Time get_clock() const;

  inline void add_cell(const std::string& species_name, const uint16_t& x, const uint16_t& y);

  List get_counts() const;

  inline List get_cells() const;

  List get_cell(const Races::Drivers::Simulation::AxisPosition& x,
                const Races::Drivers::Simulation::AxisPosition& y) const;

  List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner) const;

  List get_cells(const SEXP& first_param, const SEXP& second_param) const;

  List get_cells(const std::vector<std::string>& species_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                 const std::vector<std::string>& genotype_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  void schedule_genotype_mutation(const std::string& source, const std::string& destination,
                                  const Races::Time& time);

  void run_up_to_time(const Races::Time& time);

  void run_up_to_size(const std::string& species_name, const size_t& num_of_cells);

  void run_up_to_event(const std::string& event, const std::string& species_name,
                       const size_t& num_of_events);

  List get_firings() const;

  CharacterVector get_species_names() const;

  inline
  std::string get_directory() const;

  inline
  const std::string& get_tissue_name() const;

  IntegerVector get_tissue_size() const;

  List get_rates(const std::string& species_name) const;

  void update_rates(const std::string& species_name, const List& list);

  List choose_cell_in(const std::string& genotype_name);

  List choose_cell_in(const std::string& genotype_name,
                      const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                      const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner);

  void mutate_progeny(const Races::Drivers::Simulation::AxisPosition& x,
                      const Races::Drivers::Simulation::AxisPosition& y,
                      const std::string& mutated_genotype);

  void mutate_progeny(const List& cell_position, const std::string& mutated_genotype);

  size_t get_death_activation_level() const
  {
    return death_activation_level;
  }

  void set_death_activation_level(size_t death_activation_level)
  {
    Races::Drivers::Simulation::Simulation::death_activation_level = death_activation_level;
  }
};

Simulation::Simulation():
  Races::Drivers::Simulation::Simulation()
{}

Simulation::Simulation(const int& seed):
  Races::Drivers::Simulation::Simulation(seed)
{}

Simulation::Simulation(const std::string& output_dir):
  Races::Drivers::Simulation::Simulation(output_dir)
{}

//' @name Simulation$new
//' @title Constructs a new Simulation
//' @param output_dir The output directory (optional).
//' @param seed The seed for the pseudo-random generator (optional).
//' @examples
//' sim <- new(Simulation)
//' sim <- new(Simulation, "test")
//' sim <- new(Simulation, "test", 13)
Simulation::Simulation(const std::string& output_dir, const int& seed):
  Races::Drivers::Simulation::Simulation(output_dir, seed)
{}

//' @name Simulation$update_tissue 
//' @title Update tissue name and size
//' @param name The new name of the tissue (optional).
//' @param width The width of the new tissue.
//' @param height The height of the new tissue.
//' @examples
//' sim <- new(Simulation, "update_tissue_test")
//'
//' # set the tissue size, but not the name
//' sim$update_tissue(1200, 900)
//'
//' # set the tissue size and its name
//' sim$update_tissue("Liver", 1200, 900)
void Simulation::update_tissue(const std::string& name, const uint16_t& width, const uint16_t& height)
{
  Races::Drivers::Simulation::Simulation::set_tissue(name, {width, height});
}

void Simulation::update_tissue(const uint16_t& width, const uint16_t& height)
{
  Races::Drivers::Simulation::Simulation::set_tissue("A tissue", {width, height});
}

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
//' sim <- new(Simulation, "add_genotype_test")
//'
//' # create the two species "A+" and "A-". They both have genotype "A".
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # create the species "C" its genotype is "C".
//' sim$add_genotype(genotype = "C", growth_rate = 0.2, death_rate = 0.1)
void Simulation::add_genotype(const std::string& genotype, const List& epigenetic_rates, 
                              const List& growth_rates, const List& death_rates)
{
  using namespace Races::Drivers;

  if (!has_names(epigenetic_rates, {"+-","-+"})) {
    ::Rf_error("The second parameter must be a list specifying "
                "the epigenetic rate for \"+-\" and \"-+\"");
  }

  if (!has_names_in(growth_rates, {"+","-"})) {
    ::Rf_error("The third parameter must be a list specifying "
                "the growth rate for \"+\" and \"-\"");
  }

  if (!has_names_in(death_rates, {"+","-"})) {
    ::Rf_error("The fourth parameter must be a list specifying "
                "the death rate for \"+\" and \"-\"");
  }

  GenotypeProperties real_genotype(genotype, {{epigenetic_rates["-+"],epigenetic_rates["+-"]}});

  for (const std::string states: {"+","-"}) {
    if (growth_rates.containsElementNamed(states.c_str())) {
      real_genotype[states].set_rate(CellEventType::DUPLICATION, as<double>(growth_rates[states]));
    }
    if (death_rates.containsElementNamed(states.c_str())) {
      real_genotype[states].set_rate(CellEventType::DEATH, as<double>(death_rates[states]));
    }
  }

  Races::Drivers::Simulation::Simulation::add_genotype(real_genotype);
}

void Simulation::add_genotype(const std::string& genotype, const double& growth_rate,
                              const double& death_rate)
{
  using namespace Races::Drivers;

  GenotypeProperties real_genotype(genotype, {});

  real_genotype[""].set_rate(CellEventType::DUPLICATION, growth_rate);
  real_genotype[""].set_rate(CellEventType::DEATH, death_rate);

  Races::Drivers::Simulation::Simulation::add_genotype(real_genotype);
}

//' @name Simulation$get_species_names 
//' @title Gets the species name
//' @return The vector of the species names.
//' @examples
//' sim <- new(Simulation, "get_species_names_test")
//' sim$add_genotype("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_genotype("B", growth_rate = 0.15, death_rate = 0.05)
//'
//' # get the added species. In this case, "A" and "B"
//' sim$get_species_names()
CharacterVector Simulation::get_species_names() const
{
  CharacterVector R_species(tissue().num_of_species());

  size_t i{0};
  for (const auto& species : tissue()) {
    R_species[i] = species.get_name();
    ++i;
  }

  return R_species;
}

//' @name Simulation$add_cell 
//' @title Implant one cell in the tissue 
//' @param species The name of the new cell species.
//' @param x The position on the x axis of the cell.
//' @param y The position on the y axis of the cell.
//' @examples
//' sim <- new(Simulation, "add_cell_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # add into the tissue a cell of species "A+" in position (500,500)
//' sim$add_cell("A+", 500, 500)
void Simulation::add_cell(const std::string& species_name, const uint16_t& x, const uint16_t& y)
{
  tissue().add_cell(species_name, {x,y});
}

List Simulation::get_cells() const
{
  namespace RS = Races::Drivers::Simulation;

  std::vector<RS::AxisPosition> upper_corner = tissue().size();
  upper_corner.resize(2);

  for (auto& value : upper_corner) {
    --value;
  }

  return get_cells({0,0}, upper_corner);
}

//' @name Simulation$get_cell
//' @title Gets one of the tissue cells
//' @description This method collects some data of the aimed cell without altering 
//'      the tissue.
//' @param x The position of the aimed cell on the x axis.
//' @param y The position of the aimed cell on the y axis.
//' @return A data frame reporting "cell_id", "genotype", "epistate", "position_x", 
//'    and "position_y" of the aimed cell.
//' @examples
//' sim <- new(Simulation, "get_cell_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.3, "-" = 0.1),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)
//' sim$add_cell("A+", 500, 500)
//' sim$run_up_to_time(70)
//' 
//' # collect all the cells in the tissue
//' sim$get_cell(501, 502)
List Simulation::get_cell(const Races::Drivers::Simulation::AxisPosition& x,
                          const Races::Drivers::Simulation::AxisPosition& y) const
{
  namespace RS = Races::Drivers::Simulation;

  RS::PositionInTissue pos{x,y};

  const RS::CellInTissue& cell = tissue()({x,y});

  return wrap_a_cell(cell);
}

List Simulation::get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                           const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner) const
{
  std::set<Races::Drivers::SpeciesId> species_ids;

  for (const auto& species: tissue()) {
    species_ids.insert(species.get_id());
  }

  return get_cells(lower_corner, upper_corner, species_ids, {"+", "-"});
}

List Simulation::get_cells(const SEXP& first_param, const SEXP& second_param) const
{
  using namespace Races::Drivers::Simulation;

  if (TYPEOF(first_param)!=TYPEOF(second_param)) {
    warning(
        "The two parameters have different types: %d (%s) != %d (%s)\n",
        TYPEOF(first_param), type2name(first_param), 
        TYPEOF(second_param), type2name(second_param) 
    );
    return R_NilValue;
  }

  switch (TYPEOF(first_param)) {
      case INTSXP:
      case REALSXP:
      {
          return get_cells(as<std::vector<AxisPosition>>(first_param),
                            as<std::vector<AxisPosition>>(second_param));
      }
      case STRSXP: {
          return get_cells(as<std::vector<std::string>>(first_param),
                            as<std::vector<std::string>>(second_param));
      }
      default: {
          warning(
              "Invalid SEXPTYPE %d (%s).\n",
              TYPEOF(first_param), type2name(first_param)
          );
          return R_NilValue;
      }
  }
}

List Simulation::get_cells(const std::vector<std::string>& species_filter,
                           const std::vector<std::string>& epigenetic_filter) const
{
  namespace RS = Races::Drivers::Simulation;

  std::vector<RS::AxisPosition> upper_corner = tissue().size();
  upper_corner.resize(2);

  for (auto& value : upper_corner) {
    --value;
  }

  return get_cells({0,0}, upper_corner, species_filter, epigenetic_filter);
}

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
//' sim <- new(Simulation, "get_cells_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.02, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.3, "-" = 0.1),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)
//' sim$add_cell("A+", 500, 500)
//' sim$run_up_to_time(70)
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
List Simulation::get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                           const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                           const std::vector<std::string>& genotype_filter,
                           const std::vector<std::string>& epigenetic_filter) const
{
  std::set<std::string> genotype_set(genotype_filter.begin(), genotype_filter.end());
  std::set<std::string> epigenetic_set(epigenetic_filter.begin(), epigenetic_filter.end());

  auto species_ids = get_species_ids_from_genotype_name(tissue(), genotype_set);

  return get_cells(lower_corner, upper_corner, species_ids, epigenetic_set);
}

//' @name Simulation$get_counts
//' @title Counts the number of cells
//' @return A data frame reporting "genotype", "epistate", "counts" for each 
//'      species in the simulation.
//' @examples
//' sim <- new(Simulation, "get_counts_test")
//' sim$add_genotype("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_genotype("B", growth_rate = 0.15, death_rate = 0.05)
//' sim$schedule_genotype_mutation(src = "A", dst = "B", time = 50)
//' sim$add_cell("A", 500, 500)
//' sim$run_up_to_time(70)
//'
//' # counts the number of cells per species
//' sim$get_counts()
List Simulation::get_counts() const
{
  using namespace Races::Drivers;

  size_t num_of_rows = tissue().num_of_species();

  CharacterVector genotype_names(num_of_rows);
  CharacterVector epi_stati(num_of_rows);
  IntegerVector counts(num_of_rows);

  size_t i{0};
  for (const auto& species: tissue()) {
    genotype_names[i] = species.get_genotype_name();
    const auto& signature = species.get_methylation_signature();
    epi_stati[i] = GenotypeProperties::signature_to_string(signature);
    counts[i] = species.num_of_cells();
    ++i;
  }

  return DataFrame::create(_["genotype"]=genotype_names, _["epistate"]=epi_stati,
                            _["counts"]=counts);
}

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
//' sim <- new(Simulation, "schedule_genotype_mutation_test")
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
void Simulation::schedule_genotype_mutation(const std::string& src, const std::string& dest,
                                            const Races::Time& time)
{
  Races::Drivers::Simulation::Simulation::schedule_genotype_mutation(src,dest,time);
}

//' @name Simulation$run_up_to_time
//' @title Simulates cell evolution
//' @param time The final simulation time.
//' @examples
//' sim <- new(Simulation, "run_up_to_time_test")
//' sim$add_genotype("A", growth_rate = 0.2, death_rate = 0.1)
//' sim$add_cell("A", 500, 500)
//'
//' # simulate the tissue up to simulate timed 100
//' sim$run_up_to_time(100)
void Simulation::run_up_to_time(const Races::Time& time)
{
  Races::UI::ProgressBar bar;

  RTest<Races::Drivers::Simulation::TimeTest> ending_test{time};

  Races::Drivers::Simulation::Simulation::run(ending_test, bar);
}

//' @name Simulation$run_up_to_size
//' @title Simulates cell evolution 
//' @description This method simulates cell evolution until the number of cells in 
//'       a species reaches a specified threshold.
//' @param species The species whose number of cells is considered.
//' @param num_of_cells The threshold for the cell number.
//' @examples
//' sim <- new(Simulation, "run_up_to_size_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_cell("A+", 500, 500)
//'
//' # simulate the tissue until the species "A+" account for 100 
//' # contemporary cells
//' sim$run_up_to_size(species = "A+", num_of_cells = 100)
void Simulation::run_up_to_size(const std::string& species_name, const size_t& num_of_cells)
{
  Races::UI::ProgressBar bar;

  const auto& species_id = tissue().get_species(species_name).get_id();

  RTest<Races::Drivers::Simulation::SpeciesCountTest> ending_test{species_id, num_of_cells};

  Races::Drivers::Simulation::Simulation::run(ending_test, bar);
}

//' @name Simulation$run_up_to_event
//' @title Simulates cell evolution
//' @description This method simulates cell evolution until the number of events that 
//'         have occurred to cells of a species reaches a specified threshold.
//' @param event The considered event, i.e., "growth", "death", or "switch".
//' @param species The species whose event number is considered.
//' @param num_of_events The threshold for the event number.
//' @examples
//' sim <- new(Simulation, "run_up_to_event_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_cell("A+", 500, 500)
//'
//' # simulate the cell evolution until the number of epigenetic events from 
//' # the species "A+" is less than 100.
//' sim$run_up_to_event(event = "switch", species = "A+", num_of_events = 100)
void Simulation::run_up_to_event(const std::string& event, const std::string& species_name,
                                 const size_t& num_of_events)
{
  Races::UI::ProgressBar bar;

  if (event_names.count(event)==0) {
    handle_unknown_event(event);
  }

  namespace RS = Races::Drivers::Simulation;

  const auto& species_id = tissue().get_species(species_name).get_id();

  RTest<RS::EventCountTest> ending_test{event_names.at(event), species_id, num_of_events};

  RS::Simulation::run(ending_test, bar);
}

//' @name Simulation$get_clock 
//' @title Gets the simulated time
//' @return The time simulated by the simulation.
//' @examples
//' sim <- new(Simulation, "get_clock_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the simulated time
//' sim$get_clock()
Races::Time Simulation::get_clock() const
{
  return Races::Drivers::Simulation::Simulation::get_time();
}

//' @name Simulation$get_firings 
//' @title Gets the number of fired events 
//' @return A data frame reporting "event", "genotype", "epistate", and "fired"
//'     for each event type, genotype, and epigenetic states.
//' @examples
//' sim <- new(Simulation, "get_firings_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_cell("A+", 500, 500)
//' sim$run_up_to_event("switch", "A+", 100)
//'
//' # get the number of event fired per event and species
//' sim$get_firings()
List Simulation::get_firings() const
{
  using namespace Races::Drivers;

  size_t num_of_rows(event_names.size()*tissue().num_of_species());

  CharacterVector events(num_of_rows);
  CharacterVector genotype_names(num_of_rows);
  CharacterVector epistati(num_of_rows);
  IntegerVector firings(num_of_rows);

  const auto& t_stats = get_statistics();

  size_t i{0};
  for (const auto& species: tissue()) {
    for (const auto& [event_name, event_code]: event_names) {
      events[i] = event_name;
      genotype_names[i] = species.get_genotype_name();
      
      const auto& signature = species.get_methylation_signature();
      epistati[i] = GenotypeProperties::signature_to_string(signature);

      if (t_stats.contains_data_for(species)) {
        firings[i] = count_events(t_stats.at(species), event_code);
      } else {
        firings[i] = 0;
      }
      ++i;
    }
  }

  return DataFrame::create(_["event"]=events, _["genotype"]=genotype_names,
                            _["epistate"]=epistati, _["fired"]=firings);
}

//' @name Simulation$get_directory 
//' @title Gets the simulation directory
//' @return The directory in which the simulation is saving its progress.
//' @examples
//' sim <- new(Simulation, "test")
//'
//' # Expecting "test"
//' sim$get_directory()
std::string Simulation::get_directory() const
{ 
  return std::string(get_logger().get_directory());
}

//' @name Simulation$get_tissue_name 
//' @title Gets the tissue name
//' @return The name of the simulated tissue.
//' @examples
//' sim <- new(Simulation, "test")
//' sim$update_tissue("Liver", 1200, 900)
//'
//' # get the tissue name, i.e., expecting "Liver"
//' sim$get_tissue_name()
const std::string& Simulation::get_tissue_name() const
{
  return tissue().get_name();
}

//' @name Simulation$get_tissue_size 
//' @title Gets the size of the simulated tissue
//' @return The vector `c(x_size, y_size)` of the simulated tissue.
//' @examples
//' sim <- new(Simulation, "get_tissue_size_test")
//' sim$update_tissue("Liver", 1200, 900)
//'
//' # get the tissue size, i.e., expecting c(1200,900)
//' sim$get_tissue_size()
IntegerVector Simulation::get_tissue_size() const
{
  auto size_vect = tissue().size();

  return {size_vect[0], size_vect[1]};
}

//' @name Simulation$get_rates
//' @title Get the rates of a species
//' @param species The species whose rates are aimed.
//' @return The list of the species rates.
//' @examples
//' sim <- new(Simulation, "get_rates_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.02),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Get the rates of "A-". In this case c("growth"=0.08, "death"=0.01, "switch"=0.02) is expected
//' sim$get_rates("A-")
List Simulation::get_rates(const std::string& species_name) const
{
  auto& species = tissue().get_species(species_name);

  List rates = List::create(_("growth") = species.get_rate(event_names.at("growth")),
                            _["death"] = species.get_rate(event_names.at("death")));

  if (species.get_methylation_signature().size()>0) {
    rates.push_back(species.get_rate(event_names.at("switch")),"switch");
  }

  return rates;
}

//' @name Simulation$update_rates
//' @title Update the rates of a species
//' @param species The species whose rates must be updated.
//' @param rates The list of rates to be updated.
//' @examples
//' sim <- new(Simulation, "update_rates_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # Set the death and epigenetic switch rates of "A-" to 0
//' sim$update_rates("A-", c(switch=0, death=0))
void Simulation::update_rates(const std::string& species_name, const List& rates)
{
  using namespace Races::Drivers;
  auto& species = tissue().get_species(species_name);

  for (const auto& [event_name, event_code]: event_names) {
    if (rates.containsElementNamed(event_name.c_str())) {
      species.set_rate(event_code, as<double>(rates[event_name]));
    }
  }
}

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
//' sim <- new(Simulation, "choose_cell_in_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.1, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.15, "-" = 0.3),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_cell("A+", 500, 500)
//' sim$death_activation_level = 100
//' sim$schedule_genotype_mutation("A","B",20)
//' sim$run_up_to_size(species = "B-", num_of_cells = 50)
//'
//' # Randomly choose one cell in "B" in the tissue
//' sim$choose_cell_in(genotype = "B")
List Simulation::choose_cell_in(const std::string& genotype_name, 
                                const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                                const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner)
{
  namespace RS = Races::Drivers::Simulation;

  auto rectangle = get_rectangle(lower_corner, upper_corner);
  const auto& cell = RS::Simulation::choose_cell_in(genotype_name, rectangle);

  return wrap_a_cell(cell);
}

List Simulation::choose_cell_in(const std::string& genotype_name)
{
  namespace RS = Races::Drivers::Simulation;

  const auto& cell = RS::Simulation::choose_cell_in(genotype_name);

  return wrap_a_cell(cell);
}

void Simulation::mutate_progeny(const Races::Drivers::Simulation::AxisPosition& x,
                                const Races::Drivers::Simulation::AxisPosition& y,
                                const std::string& mutated_genotype)
{
  auto pos_in_tissue = get_position_in_tissue({x,y});

  namespace RS = Races::Drivers::Simulation;

  RS::Simulation::simulate_genotype_mutation(pos_in_tissue, mutated_genotype);
}

//' @name Simulation$mutate_progeny
//' @title Generate a mutated progeny
//' @description This method simulates both the duplication of the cell in the 
//'       specified position and the birth of one cells of a given 
//'       genotype that preserves the epigenetic status of the original cell.
//'       The mutated cell will be located in the position of its parent.
//' @param cell_position The position of the cell whose offspring will mutate.
//' @param mutated_genotype The genotype of the mutated cell.
//' @examples
//' sim <- new(Simulation, "mutate_progeny_test")
//' sim$add_genotype(genotype = "A",
//'                  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.2, "-" = 0.08),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//' sim$add_cell("A+", 500, 500)
//' sim$run_up_to_time(70)
//'
//' sim$add_genotype(genotype = "B",
//'                  epigenetic_rates = c("+-" = 0.1, "-+" = 0.01),
//'                  growth_rates = c("+" = 0.15, "-" = 0.3),
//'                  death_rates = c("+" = 0.1, "-" = 0.01))
//'
//' # duplicate the cell in position (503, 492). One of 
//' # its direct descendents will have genotype "B"
//' sim$mutate_progeny(503, 492, "B")
//'
//' # the output of `choose_cell_in` and `get_cell` can also be used
//' # as input for `mutate_progeny`
//' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
void Simulation::mutate_progeny(const List& cell_position, 
                                const std::string& mutated_genotype)
{
  namespace RS = Races::Drivers::Simulation;

  std::vector<Races::Drivers::Simulation::AxisPosition> vector_position;

  for (const std::string axis : {"x", "y"}) {
    auto field = "position_"+axis;
    if (!cell_position.containsElementNamed(field.c_str())) {
      throw std::domain_error("Missing \"" + field + "\" element from the list.");
    }
    vector_position.push_back(as<RS::AxisPosition>(cell_position[field]));
  }

  return mutate_progeny(vector_position[0], vector_position[1], mutated_genotype);
}

namespace RS = Races::Drivers::Simulation;
namespace RD = Races::Drivers;

RCPP_EXPOSED_CLASS(Simulation)
RCPP_MODULE(Drivers){
  class_<Simulation>("Simulation")
  .constructor("Create a simulation whose output file has the format \"races_<year>_<hour><minute><second>\"")
  .constructor<std::string>("Crete a simulation whose parameter is the output directory")
  .constructor<std::string, int>("Crete a simulation: the first parameter is the output directory; the second one is the random seed")

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

  // get_species_name
  .method("get_species_names", &Simulation::get_species_names,
          "Get the names of the species added to the simulation")

  // add_cell
  .method("add_cell", &Simulation::add_cell, "Place a cell in the tissue")

  // death_activation_level
  .property( "death_activation_level", &Simulation::get_death_activation_level, &Simulation::set_death_activation_level, 
          "The number of cells in a species that activate cell death" )

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
  .method("get_directory", &Simulation::get_directory, "Get the simulation directory")

  // get_tissue_name
  .method("get_tissue_name", &Simulation::get_tissue_name, "Get the simulation tissue name")

  // get_tissue_size
  .method("get_tissue_size", &Simulation::get_tissue_size, "Get the simulation tissue size")

  // get_counts
  .method("get_counts", &Simulation::get_counts, "Get the current number of cells per species")

  // get_firings
  .method("get_firings", &Simulation::get_firings,
          "Get the current number of simulated events per species")

  // get_rates
  .method("get_rates", &Simulation::get_rates, 
          "Get the rates of a species")
  
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

  // update rates
  .method("update_rates", &Simulation::update_rates, 
          "Update the rates of a species")

  // update_tissue
  .method("update_tissue", (void (Simulation::*)(const std::string&, const uint16_t&, 
                                              const uint16_t&))(&Simulation::update_tissue),
          "Update tissue name and size")
  .method("update_tissue", (void (Simulation::*)(const uint16_t&, const uint16_t&))(&Simulation::update_tissue),
          "Update tissue size");
}

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

    return static_cast<const SIMULATION_TEST*>(this)->operator()(simulation);
  }
};

const std::map<std::string, Races::Drivers::CellEventType> event_names{
  {"death",  Races::Drivers::CellEventType::DEATH},
  {"growth", Races::Drivers::CellEventType::DUPLICATION},
  {"switch", Races::Drivers::CellEventType::EPIGENETIC_EVENT},
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
    case Races::Drivers::CellEventType::EPIGENETIC_EVENT:
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

const Races::Drivers::GenotypeId& 
get_genotype_id(const Races::Drivers::Simulation::Tissue& tissue,
                const std::string& genotype_name)
{
  for (const auto& species: tissue) {
    if (species.get_genomic_name() == genotype_name) {
      return species.get_genomic_id();
    }
  }

  throw std::domain_error("Unknown genotype \""+genotype_name+"\"");
}

std::set<Races::Drivers::EpigeneticGenotypeId> 
get_epigenetic_ids_from_genotype_name(const Races::Drivers::Simulation::Tissue& tissue,
                                      const std::set<std::string>& genotype_name)
{
  std::set<Races::Drivers::EpigeneticGenotypeId > epigenetic_ids;

  for (const auto& species: tissue) {
    if (genotype_name.count(species.get_genomic_name())>0) {
      epigenetic_ids.insert(species.get_id());
    }
  }

  return epigenetic_ids;
}

size_t count_driver_mutated_cells(const Races::Drivers::Simulation::Tissue& tissue,
                                  const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                                  const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                                  const std::set<Races::Drivers::EpigeneticGenotypeId>& species_filter,
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
      if (cell_proxy.has_driver_mutations()) {
        const CellInTissue& cell = cell_proxy;

        if (species_filter.count(cell.get_epigenetic_id())>0) {

          const auto& species = tissue.get_species(cell.get_epigenetic_id());
          const auto& signature = species.get_methylation_signature();
          auto sign_string = Genotype::signature_to_string(signature);

          if (epigenetic_filter.count(sign_string)>0) {
            ++total;
          }
        }
      }
    }
  }

  return total;
}

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
                 const std::set<Races::Drivers::EpigeneticGenotypeId> &species_filter,
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
        if(cell_proxy.has_driver_mutations()) {

          const RS::CellInTissue& cell = cell_proxy;
          
          const auto& species = tissue().get_species(cell.get_epigenetic_id());
          const auto& signature = species.get_methylation_signature();
          auto sign_string = Genotype::signature_to_string(signature);
          
          if (species_filter.count(cell.get_epigenetic_id())>0
               && epigenetic_filter.count(sign_string)>0) { 

            ids[i] = cell.get_id();
            genotype_names[i] = species.get_genomic_name();

            epi_stati[i] = Genotype::signature_to_string(signature);

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

public:
  Simulation():
    Races::Drivers::Simulation::Simulation()
  {}

  Simulation(const int& seed):
    Races::Drivers::Simulation::Simulation(seed)
  {}

  Simulation(const std::string& output_dir):
    Simulation()
  {
    this->rename_log_directory(output_dir);
  }

  Simulation(const std::string& output_dir, const int& seed):
    Simulation(seed)
  {
    this->rename_log_directory(output_dir);
  }

  inline
  void set_tissue(const std::string& name, const uint16_t& width, const uint16_t& height)
  {
    static_cast<Races::Drivers::Simulation::Simulation*>(this)->set_tissue(name, {width, height});
  }

  inline
  void set_tissue(const uint16_t& width, const uint16_t& height)
  {
    static_cast<Races::Drivers::Simulation::Simulation*>(this)->set_tissue("A tissue", {width, height});
  }

  void add_species(const std::string& name, const List& epigenetic_rates, 
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

    Genotype genotype(name, {{epigenetic_rates["+-"],epigenetic_rates["-+"]}});

    for (const std::string states: {"+","-"}) {
      if (growth_rates.containsElementNamed(states.c_str())) {
        genotype[states].set_rate(CellEventType::DUPLICATION, as<double>(growth_rates[states]));
      }
      if (death_rates.containsElementNamed(states.c_str())) {
        genotype[states].set_rate(CellEventType::DEATH, as<double>(death_rates[states]));
      }
    }

    static_cast<Races::Drivers::Simulation::Simulation*>(this)->add_species(genotype);
  }

  void add_species(const std::string& name, const double& growth_rate, const double& death_rate)
  {
    using namespace Races::Drivers;

    Genotype genotype(name, {});

    genotype[""].set_rate(CellEventType::DUPLICATION, growth_rate);
    genotype[""].set_rate(CellEventType::DEATH, death_rate);

    static_cast<Races::Drivers::Simulation::Simulation*>(this)->add_species(genotype);
  }

  inline
  Races::Time get_time() const
  {
    return static_cast<const Races::Drivers::Simulation::Simulation*>(this)->get_time();
  }

  inline
  void add_cell(const std::string& epigenotype_name, const uint16_t& x, const uint16_t& y)
  {
    this->tissue().add_cell(epigenotype_name, {x,y});
  }

  List get_counts() const
  {
    using namespace Races::Drivers;

    size_t num_of_rows = tissue().num_of_species();

    CharacterVector genotype_names(num_of_rows);
    CharacterVector epi_stati(num_of_rows);
    IntegerVector counts(num_of_rows);

    size_t i{0};
    for (const auto& species: tissue()) {
      genotype_names[i] = species.get_genomic_name();
      const auto& signature = species.get_methylation_signature();
      epi_stati[i] = Genotype::signature_to_string(signature);
      counts[i] = species.num_of_cells();
      ++i;
    }

    return DataFrame::create(_["genotype"]=genotype_names, _["epistate"]=epi_stati,
                             _["counts"]=counts);
  }

  inline List get_cells() const
  {
    namespace RS = Races::Drivers::Simulation;
    
    std::vector<RS::AxisPosition> upper_corner = tissue().size();
    upper_corner.resize(2);

    for (auto& value : upper_corner) {
      --value;
    }

    return get_cells({0,0}, upper_corner);
  }

  List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner) const
  {
    std::set<Races::Drivers::EpigeneticGenotypeId> species_ids;

    for (const auto& species: tissue()) {
      species_ids.insert(species.get_id());
    }

    return get_cells(lower_corner, upper_corner, species_ids, {"+", "-"});
  }

  List get_cells(const SEXP& first_param, const SEXP& second_param) const
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

  List get_cells(const std::vector<std::string>& species_filter,
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

  List get_cells(const std::vector<Races::Drivers::Simulation::AxisPosition>& lower_corner, 
                 const std::vector<Races::Drivers::Simulation::AxisPosition>& upper_corner,
                 const std::vector<std::string>& genotype_filter,
                 const std::vector<std::string>& epigenetic_filter) const
  {
    std::set<std::string> genotype_set(genotype_filter.begin(), genotype_filter.end());
    std::set<std::string> epigenetic_set(epigenetic_filter.begin(), epigenetic_filter.end());

    auto species_ids = get_epigenetic_ids_from_genotype_name(tissue(), genotype_set);

    return get_cells(lower_corner, upper_corner, species_ids, epigenetic_set);
  }

  void add_driver_mutation(const std::string& source, const std::string& destination,
                           const Races::Time& time)
  {
    namespace RS = Races::Drivers::Simulation;

    static_cast<RS::Simulation*>(this)->add_driver_mutation(source,destination,time);
  }

  void run_up_to_time(const Races::Time& time)
  {
    Races::UI::ProgressBar bar;

    RTest<Races::Drivers::Simulation::TimeTest> ending_test{time};

    static_cast<Races::Drivers::Simulation::Simulation*>(this)->run(ending_test, bar);
  }

  void run_up_to_size(const std::string& species_name, const size_t& num_of_cells)
  {
    Races::UI::ProgressBar bar;

    const auto& species_id = tissue().get_species(species_name).get_id();

    RTest<Races::Drivers::Simulation::SpeciesCountTest> ending_test{species_id, num_of_cells};

    static_cast<Races::Drivers::Simulation::Simulation*>(this)->run(ending_test, bar);
  }

  void run_up_to_event(const std::string& event, const std::string& genomic_name,
                       const std::string& epistate, const size_t& threshold)
  {
    Races::UI::ProgressBar bar;

    if (event_names.count(event)==0) {
      handle_unknown_event(event);
    }

    namespace RS = Races::Drivers::Simulation;

    const auto& species_id = tissue().get_species(genomic_name+epistate).get_id();

    RTest<RS::EventCountTest> ending_test{event_names.at(event), species_id, threshold};

    static_cast<RS::Simulation*>(this)->run(ending_test, bar);
  }

  List get_firings() const
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
        genotype_names[i] = species.get_genomic_name();
        
        const auto& signature = species.get_methylation_signature();
        epistati[i] = Genotype::signature_to_string(signature);

        if (t_stats.contains_data_for(species)) {
          firings[i] = count_events(t_stats.at(species), event_code);
        } else {
          firings[i] = 0;
        }
        ++i;
      }
    }

    return DataFrame::create(_["events"]=events, _["genotype"]=genotype_names,
                             _["epistate"]=epistati, _["firings"]=firings);
  }

  List get_species_names() const
  {
    auto genotypes = this->tissue().get_genotypes();

    List R_genotypes(genotypes.size());

    for (size_t i=0; i<genotypes.size(); ++i) {
      R_genotypes[i] = genotypes[i].get_epigenetic_name();
    }

    return R_genotypes;
  }

  inline
  std::string get_directory() const
  { 
    return std::string(get_logger().get_directory());
  }

  inline
  const std::string& get_tissue_name() const
  {
    return tissue().get_name();
  }

  IntegerVector get_tissue_size() const
  {
    auto size_vect = tissue().size();

    return {size_vect[0], size_vect[1]};
  }
};

namespace RS = Races::Drivers::Simulation;
namespace RD = Races::Drivers;

RCPP_EXPOSED_CLASS(Simulation)
RCPP_MODULE(Drivers){
  class_<Simulation>("Simulation")
  .constructor()
  .constructor<std::string>()
  .constructor<std::string, int>()
  
  // set_tissue
  .method("set_tissue", (void (Simulation::*)(const std::string&, const uint16_t&, 
                                              const uint16_t&))(&Simulation::set_tissue),
          "Set tissue name and size")
  .method("set_tissue", (void (Simulation::*)(const uint16_t&, const uint16_t&))(&Simulation::set_tissue),
          "Set tissue size")

  // add_species
  .method("add_species", (void (Simulation::*)(const std::string&, const List&, const List&,
                                               const List&))(&Simulation::add_species),
          "Add a new species with epigenetic status")
  .method("add_species", (void (Simulation::*)(const std::string&, const double&,
                                               const double&))(&Simulation::add_species),
          "Add a new species")

  // add_timed_mutation
  .method("add_timed_mutation", &Simulation::add_driver_mutation,
          "Add a timed mutation between two different species")

  // get_species_name
  .method("get_species_names", &Simulation::get_species_names,
          "Get the names of the species added to the simulation")

  // add_cell
  .method("add_cell", &Simulation::add_cell, "Place a cell in the tissue")

  // get_clock
  .method("get_clock", &Simulation::get_time, "Get the current simulation time")

  // get_cells
  .method("get_cells", (List (Simulation::*)(const std::vector<RS::AxisPosition>&, 
                                             const std::vector<RS::AxisPosition>&,
                                             const std::vector<std::string>&,
                                             const std::vector<std::string>&) const)(&Simulation::get_cells),
          "Get cell data from the simulated tissue")
  .method("get_cells", (List (Simulation::*)(const SEXP&, const SEXP&) const)(&Simulation::get_cells),
          "Get cell data from the simulated tissue")
  .method("get_cells", (List (Simulation::*)() const)(&Simulation::get_cells),
          "Get cell data from the simulated tissue")

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

  // run_up_to_time
  .method("run_up_to_time", &Simulation::run_up_to_time, 
          "Simulate the system up to the specified simulation time")

  // run_up_to_event
  .method("run_up_to_event", &Simulation::run_up_to_event, 
          "Simulate the system up to the specified number of events")

  // run_up_to_size
  .method("run_up_to_size", &Simulation::run_up_to_size, 
          "Simulate the system up to the specified number of cells in the species");
}

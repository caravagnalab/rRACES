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

#include <algorithm>
#include <filesystem>

#include <ending_conditions.hpp>

#include "simulation.hpp"

#include "utility.hpp"


template<typename SIMULATION_TEST>
struct RTest : public SIMULATION_TEST
{
  size_t counter;

  template<typename ...Args>
  explicit RTest(Args...args):
      SIMULATION_TEST(args...), counter{0}
  {}

  bool operator()(const RACES::Mutants::Evolutions::Simulation& simulation)
  {
    if (++counter >= 10000) {
      counter = 0;
      try {
        Rcpp::checkUserInterrupt();
      } catch (...) {
        return true;
      }
    }

    using namespace RACES::Mutants::Evolutions;

    return SIMULATION_TEST::operator()(simulation);
  }
};

const std::map<std::string, RACES::Mutants::CellEventType> event_names{
  {"death",  RACES::Mutants::CellEventType::DEATH},
  {"growth", RACES::Mutants::CellEventType::DUPLICATION},
  {"switch", RACES::Mutants::CellEventType::EPIGENETIC_SWITCH},
};

size_t count_events(const RACES::Mutants::Evolutions::SpeciesStatistics& statistics,
                    const RACES::Mutants::CellEventType& event)
{
  switch(event)
  {
    case RACES::Mutants::CellEventType::DEATH:
      return statistics.killed_cells;
    case RACES::Mutants::CellEventType::DUPLICATION:
      return statistics.num_duplications;
    case RACES::Mutants::CellEventType::EPIGENETIC_SWITCH:
      return statistics.num_of_epigenetic_events();
    default:
      ::Rf_error("get_counts: unsupported event");
  }
}

inline std::string get_signature_string(const RACES::Mutants::Evolutions::Species& species)
{
  const auto& signature = species.get_methylation_signature();
  return RACES::Mutants::MutantProperties::signature_to_string(signature);
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
      oss << " ";
    }

    if ((++i)==event_names.size()) {
      oss << "and ";
    }

    oss << "\"" << name << "\"";
  }

  oss << ".";

  throw std::domain_error(oss.str());
}

std::set<RACES::Mutants::SpeciesId>
get_species_ids_from_mutant_name(const RACES::Mutants::Evolutions::Tissue& tissue,
                                 const std::set<std::string>& mutant_name)
{
  std::set<RACES::Mutants::SpeciesId > species_ids;

  for (const auto& species: tissue) {
    if (mutant_name.count(species.get_mutant_name())>0) {
      species_ids.insert(species.get_id());
    }
  }

  return species_ids;
}

RACES::Mutants::Evolutions::PositionInTissue
get_position_in_tissue(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& position)
{
  if (position.size()==2) {
    return {position[0], position[1]};
  }

  ::Rf_error("rRACES supports only 2 dimensional space so far");
}

RACES::Mutants::RectangleSet
get_rectangle(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
              const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner)
{
  auto l_position = get_position_in_tissue(lower_corner);
  auto u_position = get_position_in_tissue(upper_corner);

  return {l_position, u_position};
}

size_t count_driver_mutated_cells(const RACES::Mutants::Evolutions::Tissue& tissue,
                                  const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                                  const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                                  const std::set<RACES::Mutants::SpeciesId>& species_filter,
                                  const std::set<std::string>& epigenetic_filter)
{
  using namespace RACES::Mutants;
  using namespace RACES::Mutants::Evolutions;

  if (lower_corner.size() != upper_corner.size()) {
    ::Rf_error("lower_corner and upper_corner must have the same size");
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
          auto sign_string = get_signature_string(species);

          if (epigenetic_filter.count(sign_string)>0) {
            ++total;
          }
        }
      }
    }
  }

  return total;
}

std::vector<RACES::Mutants::Evolutions::Direction> SpatialSimulation::get_possible_directions()
{
  namespace RS = RACES::Mutants::Evolutions;

  std::vector<RS::Direction> directions;
  for (const auto &x_move : {RS::Direction::X_UP, RS::Direction::X_DOWN, RS::Direction::X_NULL}) {
      for (const auto &y_move : {RS::Direction::Y_UP, RS::Direction::Y_DOWN, RS::Direction::Y_NULL}) {
          directions.push_back(x_move|y_move);
      }
  }

  // remove null move
  directions.pop_back();

  return directions;
}


PlainChooser::PlainChooser(const std::shared_ptr<RACES::Mutants::Evolutions::Simulation>& sim_ptr,
                           const std::string& mutant_name):
  sim_ptr(sim_ptr), mutant_name(mutant_name)
{}


RectangularChooser::RectangularChooser(
          const std::shared_ptr<RACES::Mutants::Evolutions::Simulation>& sim_ptr,
          const std::string& mutant_name,
          const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
          const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner):
  PlainChooser(sim_ptr, mutant_name), rectangle(get_rectangle(lower_corner, upper_corner))
{}

bool SpatialSimulation::has_names(const Rcpp::List& list, std::vector<std::string> aimed_names)
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

bool SpatialSimulation::has_names_in(const Rcpp::List& list, std::set<std::string> aimed_names)
{
  using namespace Rcpp;

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

Rcpp::List
SpatialSimulation::get_cells(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                      const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                      const std::set<RACES::Mutants::SpeciesId> &species_filter,
                      const std::set<std::string> &epigenetic_filter) const
{
  using namespace Rcpp;
  using namespace RACES::Mutants;

  namespace RS = RACES::Mutants::Evolutions;

  if (lower_corner.size() != 2) {
    ::Rf_error("The lower corner must be a vector having size 2");
  }

  if (upper_corner.size() != 2) {
    ::Rf_error("The upper corner must be a vector having size 2");
  }

  size_t num_of_rows = count_driver_mutated_cells(sim_ptr->tissue(), lower_corner, upper_corner,
                                                  species_filter, epigenetic_filter);

  IntegerVector ids(num_of_rows);
  CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);
  IntegerVector x_pos(num_of_rows), y_pos(num_of_rows);

  size_t i{0};
  for (auto x=lower_corner[0]; x<=upper_corner[0]; ++x) {
    for (auto y=lower_corner[1]; y<=upper_corner[1]; ++y) {
      auto cell_proxy = sim_ptr->tissue()({x,y});
      if(!cell_proxy.is_wild_type()) {

        const RS::CellInTissue& cell = cell_proxy;

        const auto& species = sim_ptr->tissue().get_species(cell.get_species_id());
        const auto sign_string = get_signature_string(species);

        if (species_filter.count(cell.get_species_id())>0
              && epigenetic_filter.count(sign_string)>0) {

          ids[i] = cell.get_id();
          mutant_names[i] = species.get_mutant_name();

          epi_states[i] = sign_string;

          x_pos[i] = x;
          y_pos[i] = y;

          ++i;
        }
      }
    }
  }

  return DataFrame::create(_["cell_id"]=ids, _["mutant"]=mutant_names,
                            _["epistate"]=epi_states, _["position_x"]=x_pos,
                            _["position_y"]=y_pos);
}

Rcpp::List SpatialSimulation::wrap_a_cell(const RACES::Mutants::Evolutions::CellInTissue& cell) const
{
  using namespace Rcpp;
  using namespace RACES::Mutants;

  const auto& species = sim_ptr->tissue().get_species(cell.get_species_id());

  const auto& mutant_name = sim_ptr->find_mutant_name(species.get_mutant_id());

  auto epistate = MutantProperties::signature_to_string(species.get_methylation_signature());

  return DataFrame::create(_["cell_id"]=cell.get_id(), _["mutant"]=mutant_name,
                            _["epistate"]=epistate, _["position_x"]=cell.x,
                            _["position_y"]=cell.y);
}

SpatialSimulation SpatialSimulation::load(const std::string& directory_name)
{
  using namespace RACES::Mutants::Evolutions;

  SpatialSimulation simulation;

  simulation.save_snapshots = true;
  simulation.name = directory_name;

  if (!std::filesystem::exists(directory_name)) {
    throw std::domain_error("The directory \"" + directory_name + "\" does not exist.");
  }

  if (!std::filesystem::is_directory(directory_name)) {
    throw std::domain_error("\"" + directory_name + "\" is not a directory.");
  }

  auto snapshot_path = BinaryLogger::find_last_snapshot_in(directory_name);

  RACES::Archive::Binary::In archive(snapshot_path);

  try {
    archive & *(simulation.sim_ptr);
  } catch (RACES::Archive::WrongFileFormatDescr& ex) {
    raise_error(ex, "spatial simulation");
  } catch (RACES::Archive::WrongFileFormatVersion& ex) {
    raise_error(ex, "spatial simulation");
  }

  auto ruh_path = std::filesystem::path(directory_name)/
                            get_rates_update_history_file_name();

  if (std::filesystem::exists(ruh_path)) {
    RACES::Archive::Binary::In ruh_archive(ruh_path);

    ruh_archive & simulation.rate_update_history;
  } else {
    Rcpp::warning("The rates update history file is missing.");
  }
  return simulation;
}

std::string get_time_string()
{
    std::time_t time;
    std::tm* info;
    char buffer[81];

    std::time(&time);
    info = std::localtime(&time);

    std::strftime(buffer,80,"%Y%m%d-%H%M%S",info);

    return buffer;
}

inline std::string
get_default_name()
{
  return "races_"+get_time_string();
}

void SpatialSimulation::init(const SEXP& sexp)
{
  using namespace Rcpp;

  namespace RS = RACES::Mutants::Evolutions;

  switch (TYPEOF(sexp)) {
    case INTSXP:
    case REALSXP:
    {
      int seed = as<int>(sexp);
      name = get_default_name();

      if (save_snapshots) {
        sim_ptr = std::make_shared<RS::Simulation>(name, seed);
      } else {
        sim_ptr = std::make_shared<RS::Simulation>(get_tmp_dir_path(), seed);
      }
      break;
    }
    case STRSXP: {
      name = as<std::string>(sexp);

      auto seed = get_random_seed<int>(R_NilValue);
      if (save_snapshots) {
        sim_ptr = std::make_shared<RS::Simulation>(name, seed);
      } else {
        sim_ptr = std::make_shared<RS::Simulation>(get_tmp_dir_path(), seed);
      }
      break;
    }
    default: {
      std::ostringstream oss;

      oss << "Invalid type for the first parameter: "
          << type2name(sexp);

      throw std::domain_error(oss.str());
    }
  }
}

SpatialSimulation::SpatialSimulation():
  sim_ptr(std::make_shared<RACES::Mutants::Evolutions::Simulation>(get_tmp_dir_path(),
                                                                   get_random_seed<int>(R_NilValue))),
  name(get_default_name()), save_snapshots(false)
{}

SpatialSimulation::SpatialSimulation(const SEXP& sexp):
  save_snapshots(false)
{
  using namespace Rcpp;

  if (TYPEOF(sexp) == LGLSXP) {
    save_snapshots = as<bool>(sexp);
    name = get_default_name();

    auto seed = get_random_seed<int>(R_NilValue);
    if (save_snapshots) {
      sim_ptr = std::make_shared<RACES::Mutants::Evolutions::Simulation>(name, seed);
    } else {
      sim_ptr = std::make_shared<RACES::Mutants::Evolutions::Simulation>(get_tmp_dir_path(), seed);
    }

    return;
  }

  init(sexp);
}

SpatialSimulation::SpatialSimulation(const SEXP& first_param, const SEXP& second_param):
  save_snapshots(false)
{
  using namespace Rcpp;

  if (TYPEOF(second_param) == LGLSXP) {
    save_snapshots = as<bool>(second_param);

    init(first_param);

    return;
  }

  if (TYPEOF(first_param) != STRSXP) {
    std::ostringstream oss;

    oss << "Invalid type for the parameter 1: "
        << type2name(first_param)
        << ". If the last parameter is not a Boolean value (save on disk"
        << " parameter), it must be a string (the name of the simulation).";

    throw std::domain_error(oss.str());
  }

  if (TYPEOF(second_param) != INTSXP && TYPEOF(second_param) != REALSXP) {
    std::ostringstream oss;

    oss << "Invalid type for the parameter 2: "
        << type2name(second_param)
        << ". If the last parameter is not a Boolean value (save on disk"
        << " parameter), it must be an integer value (the random seed).";

    throw std::domain_error(oss.str());
  }

  name = as<std::string>(first_param);
  int seed = as<int>(second_param);

  sim_ptr = std::make_shared<RACES::Mutants::Evolutions::Simulation>(name, seed);
}

SpatialSimulation::SpatialSimulation(const std::string& simulation_name,
                                            const SEXP& seed,
                                            const bool& save_snapshots):
  name(simulation_name), save_snapshots(save_snapshots)
{
  int c_seed = get_random_seed<int>(seed);

  if (save_snapshots) {
    sim_ptr = std::make_shared<RACES::Mutants::Evolutions::Simulation>(simulation_name, c_seed);
  } else {
    sim_ptr = std::make_shared<RACES::Mutants::Evolutions::Simulation>(get_tmp_dir_path(), c_seed);
  }
}

SpatialSimulation::SpatialSimulation(const std::string& simulation_name, const int& seed,
                                     const bool& save_snapshots):
  name(simulation_name), save_snapshots(save_snapshots)
{
  if (save_snapshots) {
    sim_ptr = std::make_shared<RACES::Mutants::Evolutions::Simulation>(simulation_name, seed);
  } else {
    sim_ptr = std::make_shared<RACES::Mutants::Evolutions::Simulation>(get_tmp_dir_path(), seed);
  }
}

std::string get_string(const SEXP& parameter, const std::string parameter_name)
{
  using namespace Rcpp;

  if (TYPEOF(parameter) != STRSXP) {
    throw std::domain_error("The parameter \"" + parameter_name
                            + "\" must be a string.");
  }

  return as<std::string>(parameter); 
}

bool get_bool(const SEXP& parameter, const std::string parameter_name)
{
  using namespace Rcpp;

  if (TYPEOF(parameter) != LGLSXP) {
    throw std::domain_error("The parameter \"" + parameter_name
                           + "\" must be a Boolean value.");
  }

  return as<int>(parameter); 
}

size_t get_size(const SEXP& parameter, const std::string parameter_name)
{
  using namespace Rcpp;

  if (TYPEOF(parameter) != INTSXP && TYPEOF(parameter) != REALSXP) {
    throw std::domain_error("The parameter \"" + parameter_name
                            + "\" must be a positive integer.");
  }

  auto c_value = as<long>(parameter); 

  if (c_value <= 0) {
    throw std::domain_error("The parameter \"" + parameter_name
                            + "\" must be a positive integer.");
  }

  return static_cast<size_t>(c_value); 
}

SpatialSimulation
SpatialSimulation::build_simulation(const SEXP& simulation_name, const SEXP& width,
                                    const SEXP& height, const SEXP& save_snapshots,
                                    const SEXP& seed)
{
  std::string c_name;
  if (TYPEOF(simulation_name) == NILSXP) {
    c_name = to_string(get_tmp_dir_path());
  } else {
    c_name = get_string(simulation_name, "name");
  }
  auto c_width = get_size(width, "width");
  auto c_height = get_size(height, "height");
  auto c_save = get_bool(save_snapshots, "save_snapshots");
  auto c_seed = get_random_seed<int>(seed);

  SpatialSimulation sim(c_name, c_seed, c_save);

  sim.update_tissue(c_width, c_height);

  return sim;
}

SpatialSimulation::~SpatialSimulation()
{
  if (sim_ptr.use_count()==1 && !save_snapshots) {
    auto dir = sim_ptr->get_logger().get_directory();

    sim_ptr = nullptr;

    std::filesystem::remove_all(dir);
  }
}

void SpatialSimulation::add_mutant_rate_history(const RACES::Mutants::MutantProperties& mutant_propeties)
{
  auto& timed_update = rate_update_history[sim_ptr->get_time()];
  for (const auto& species : mutant_propeties.get_species()) {
    auto& species_update = timed_update[species.get_id()];

    const auto event_rates = species.get_rates();
    for (const auto& [event_name, event_code]: event_names) {
        auto found = event_rates.find(event_code);

        if (found != event_rates.end()) {
            species_update[event_name] = found->second;
        }
    }
  }

  auto sim_path = sim_ptr->get_logger().get_directory();
  if (!std::filesystem::exists(sim_path)) {
    std::filesystem::create_directory(sim_path);
  }

  RACES::Archive::Binary::Out ruh_archive(sim_path/get_rates_update_history_file_name());

  ruh_archive & rate_update_history;
}

void SpatialSimulation::add_mutant(const std::string& mutant_name, const Rcpp::List& epigenetic_rates,
                            const Rcpp::List& growth_rates, const Rcpp::List& death_rates)
{
  using namespace Rcpp;
  using namespace RACES::Mutants;

  if (mutant_name.find('+')!=std::string::npos
      || mutant_name.find('-')!=std::string::npos) {
    ::Rf_error("Mutant name cannot contains a '-' or a '+'.");
  }

  if (mutant_name.find('.')!=std::string::npos) {
    ::Rf_error("Mutant name cannot contains a '.'.");
  }

  if (mutant_name == "Wild-type") {
    ::Rf_error("\"Wild-type\" is a reserved mutant name.");
  }

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

  MutantProperties real_mutant(mutant_name, {{epigenetic_rates["-+"],epigenetic_rates["+-"]}});

  for (const std::string states: {"+","-"}) {
    if (growth_rates.containsElementNamed(states.c_str())) {
      real_mutant[states].set_rate(CellEventType::DUPLICATION, as<double>(growth_rates[states]));
    }
    if (death_rates.containsElementNamed(states.c_str())) {
      real_mutant[states].set_rate(CellEventType::DEATH, as<double>(death_rates[states]));
    }
  }

  sim_ptr->add_mutant(real_mutant);

  add_mutant_rate_history(real_mutant);
}

void SpatialSimulation::add_mutant(const std::string& mutant_name, const double& growth_rate,
                            const double& death_rate)
{
  using namespace RACES::Mutants;

  if (mutant_name == "Wild-type") {
    ::Rf_error("\"Wild-type\" is a reserved mutant name.");
  }

  MutantProperties real_mutant(mutant_name, {});

  real_mutant[""].set_rate(CellEventType::DUPLICATION, growth_rate);
  real_mutant[""].set_rate(CellEventType::DEATH, death_rate);

  sim_ptr->add_mutant(real_mutant);

  add_mutant_rate_history(real_mutant);
}

Rcpp::List SpatialSimulation::get_species() const
{
  using namespace Rcpp;
  size_t num_of_rows = sim_ptr->tissue().num_of_species();

  CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);
  NumericVector switch_rates(num_of_rows), duplication_rates(num_of_rows),
                death_rates(num_of_rows);

  using namespace RACES::Mutants;

  size_t i{0};
  for (const auto& species: sim_ptr->tissue()) {
    mutant_names[i] = species.get_mutant_name();
    duplication_rates[i] = species.get_rate(CellEventType::DUPLICATION);
    death_rates[i] = species.get_rate(CellEventType::DEATH);
    epi_states[i] = get_signature_string(species);

    const auto& species_switch_rates = species.get_epigenetic_switch_rates();
    switch(species_switch_rates.size()) {
      case 0:
        switch_rates[i] = NA_REAL;
        break;
      case 1:
        switch_rates[i] = species_switch_rates.begin()->second;
        break;
      default:
        ::Rf_error("rRACES does not support multiple promoters");
    }

    ++i;
  }

  return DataFrame::create(_["mutant"]=mutant_names, _["epistate"]=epi_states,
                            _["growth_rate"]=duplication_rates,
                            _["death_rate"]=death_rates,
                            _["switch_rate"]=switch_rates);
}

Rcpp::List SpatialSimulation::get_rates_update_history() const
{
  using namespace Rcpp;
  using namespace RACES::Mutants;

  CharacterVector mutant_names, epi_states, event_names;
  NumericVector rates, times;

  const auto& tissue = sim_ptr->tissue();
  for (const auto& [time, species_rate_updates] : rate_update_history) {
    for (const auto& [species_id, event_rate_updates] : species_rate_updates) {
      const auto& species = tissue.get_species(species_id);
      const std::string mutant_name = sim_ptr->find_mutant_name(species.get_mutant_id());
      const auto& m_signature = species.get_methylation_signature();
      const auto epi_state = MutantProperties::signature_to_string(m_signature);
      for (const auto& [event_name, rate] : event_rate_updates) {
        times.push_back(time);
        mutant_names.push_back(mutant_name.c_str());
        epi_states.push_back(epi_state.c_str());
        event_names.push_back(event_name.c_str());
        rates.push_back(rate);
      }
    }
  }

  return DataFrame::create(_["time"]=times, _["mutant"]=mutant_names,
                           _["epistate"]=epi_states, _["event"]=event_names,
                           _["rate"]=rates);
}

void SpatialSimulation::place_cell(const std::string& species_name,
                            const RACES::Mutants::Evolutions::AxisPosition& x,
                            const RACES::Mutants::Evolutions::AxisPosition& y)
{
  if (sim_ptr->tissue().num_of_mutated_cells()>0) {
    Rcpp::warning("Warning: the tissue already contains a cell.");
  }

  const auto& species = sim_ptr->tissue().get_species(species_name);

  sim_ptr->place_cell(species.get_id(), {x,y});
}

Rcpp::List SpatialSimulation::get_cells() const
{
  namespace RS = RACES::Mutants::Evolutions;

  std::vector<RS::AxisPosition> upper_corner = sim_ptr->tissue().size();
  upper_corner.resize(2);

  for (auto& value : upper_corner) {
    --value;
  }

  return get_cells({0,0}, upper_corner);
}

Rcpp::List SpatialSimulation::get_cell(const RACES::Mutants::Evolutions::AxisPosition& x,
                          const RACES::Mutants::Evolutions::AxisPosition& y) const
{
  namespace RS = RACES::Mutants::Evolutions;

  const RS::CellInTissue& cell = sim_ptr->tissue()({x,y});

  return wrap_a_cell(cell);
}

Rcpp::List SpatialSimulation::get_cells(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                           const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner) const
{
  std::set<RACES::Mutants::SpeciesId> species_ids;

  for (const auto& species: sim_ptr->tissue()) {
    species_ids.insert(species.get_id());
  }

  return get_cells(lower_corner, upper_corner, species_ids, {"+", "-", ""});
}

Rcpp::List SpatialSimulation::get_cells(const SEXP& first_param, const SEXP& second_param) const
{
  using namespace Rcpp;
  using namespace RACES::Mutants::Evolutions;

  if (TYPEOF(first_param)!=TYPEOF(second_param)) {
    std::ostringstream oss;

    oss << "The two parameters have different types: " << type2name(first_param)
        << " != " << type2name(second_param);

    throw std::domain_error(oss.str());
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
        std::ostringstream oss;

        oss << "Invalid parameter type " << type2name(first_param);

        throw std::domain_error(oss.str());
      }
  }
}

Rcpp::List SpatialSimulation::get_cells(const std::vector<std::string>& species_filter,
                           const std::vector<std::string>& epigenetic_filter) const
{
  namespace RS = RACES::Mutants::Evolutions;

  std::vector<RS::AxisPosition> upper_corner = sim_ptr->tissue().size();
  upper_corner.resize(2);

  for (auto& value : upper_corner) {
    --value;
  }

  return get_cells({0,0}, upper_corner, species_filter, epigenetic_filter);
}

Rcpp::List SpatialSimulation::get_cells(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                           const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                           const std::vector<std::string>& mutant_filter,
                           const std::vector<std::string>& epigenetic_filter) const
{
  std::set<std::string> mutant_set(mutant_filter.begin(), mutant_filter.end());
  std::set<std::string> epigenetic_set(epigenetic_filter.begin(), epigenetic_filter.end());

  auto species_ids = get_species_ids_from_mutant_name(sim_ptr->tissue(), mutant_set);

  return get_cells(lower_corner, upper_corner, species_ids, epigenetic_set);
}

Rcpp::List SpatialSimulation::get_counts() const
{
  using namespace Rcpp;
  using namespace RACES::Mutants;

  size_t num_of_rows = sim_ptr->tissue().num_of_species();

  CharacterVector mutant_names(num_of_rows);
  CharacterVector epi_states(num_of_rows);
  IntegerVector counts(num_of_rows);

  size_t i{0};
  for (const auto& species: sim_ptr->tissue()) {
    mutant_names[i] = species.get_mutant_name();
    epi_states[i] = get_signature_string(species);
    counts[i] = species.num_of_cells();
    ++i;
  }

  return DataFrame::create(_["mutant"]=mutant_names, _["epistate"]=epi_states,
                            _["counts"]=counts);
}

std::map<RACES::Mutants::SpeciesId, std::string>
get_species_id2name(const RACES::Mutants::Evolutions::Tissue& tissue)
{
  std::map<RACES::Mutants::SpeciesId, std::string> id2name;
  for (const auto& species : tissue) {
    id2name[species.get_id()] = species.get_name();
  }

  return id2name;
}

Rcpp::List SpatialSimulation::get_added_cells() const
{
  using namespace Rcpp;
  using namespace RACES::Mutants;

  namespace RS = RACES::Mutants::Evolutions;

  size_t num_of_rows = sim_ptr->get_added_cells().size();

  CharacterVector mutant_names(num_of_rows),  epi_states(num_of_rows);
  IntegerVector position_x(num_of_rows), position_y(num_of_rows);
  NumericVector time(num_of_rows);

  size_t i{0};
  for (const auto& added_cell: sim_ptr->get_added_cells()) {
    const auto& species = sim_ptr->tissue().get_species(added_cell.species_id);
    mutant_names[i] = sim_ptr->find_mutant_name(species.get_mutant_id());
    epi_states[i] = get_signature_string(species);
    position_x[i] = added_cell.x;
    position_y[i] = added_cell.y;
    time[i] = added_cell.time;
    ++i;
  }

  return DataFrame::create(_["mutant"]=mutant_names, _["epistate"]=epi_states,
                           _["position_x"]=position_x,  _["position_y"]=position_y,
                           _["time"] = time);
}

// sorting LineageEdge by time
struct TimedLineageEdge : public RACES::Mutants::Evolutions::LineageEdge
{
  RACES::Time time;

  TimedLineageEdge():
    RACES::Mutants::Evolutions::LineageEdge(), time(0)
  {}

  TimedLineageEdge(const RACES::Mutants::Evolutions::LineageEdge& edge, const RACES::Time& time):
    RACES::Mutants::Evolutions::LineageEdge(edge), time(time)
  {}
};

struct TimedLineageEdgeCmp
{
  bool operator()(const TimedLineageEdge& a, const TimedLineageEdge& b)
  {
    return (a.time<b.time
            || (a.time==b.time && (a.get_ancestor()<b.get_ancestor()))
            || (a.time==b.time && (a.get_ancestor()==b.get_ancestor())
                && (a.get_progeny()<b.get_progeny())));
  }
};

std::vector<TimedLineageEdge> sorted_timed_edges(const RACES::Mutants::Evolutions::Simulation& simulation)
{
  const auto& lineage_graph = simulation.get_lineage_graph();
  const size_t num_of_edges = lineage_graph.num_of_edges();

  std::vector<TimedLineageEdge> timed_edges;

  timed_edges.reserve(num_of_edges);

  for (const auto& [edge, edge_time] : lineage_graph) {
    timed_edges.push_back({edge, edge_time});
  }

  TimedLineageEdgeCmp cmp;
  sort(timed_edges.begin(), timed_edges.end(), cmp);

  return timed_edges;
}

Rcpp::List SpatialSimulation::get_lineage_graph() const
{
  using namespace Rcpp;
  const auto species_id2name = get_species_id2name(sim_ptr->tissue());

  const auto timed_edges = sorted_timed_edges(*sim_ptr);

  CharacterVector ancestors(timed_edges.size()), progeny(timed_edges.size());
  NumericVector first_cross(timed_edges.size());

  size_t i{0};
  for (const auto& timed_edge : timed_edges) {
    ancestors[i] = (timed_edge.get_ancestor() != WILD_TYPE_SPECIES ?
                    species_id2name.at(timed_edge.get_ancestor()):
                    "Wild-type");

    progeny[i] = (timed_edge.get_progeny() != WILD_TYPE_SPECIES ?
                  species_id2name.at(timed_edge.get_progeny()):
                  "Wild-type");
    first_cross[i] = timed_edge.time;

    ++i;
  }

  return DataFrame::create(_["ancestor"]=ancestors, _["progeny"]=progeny,
                            _["first_cross"]=first_cross);
}

inline void validate_non_empty_tissue(const RACES::Mutants::Evolutions::Tissue& tissue)
{
  if (tissue.num_of_cells()==0) {
    ::Rf_error("The tissue does not contain any cell.");
  }
}

void SpatialSimulation::run_up_to_time(const RACES::Time& time)
{
  run_up_to_time(time, false);
}

void SpatialSimulation::run_up_to_time(const RACES::Time& time, const bool quiet)
{
  validate_non_empty_tissue(sim_ptr->tissue());

  RTest<RACES::Mutants::Evolutions::TimeTest> ending_test{time};

  RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

  sim_ptr->run(ending_test, progress_bar);
}

void SpatialSimulation::run_up_to_size(const std::string& species_name, const size_t& num_of_cells)
{
    run_up_to_size(species_name, num_of_cells, false);
}

void SpatialSimulation::run_up_to_size(const std::string& species_name, const size_t& num_of_cells, const bool quiet)
{
  validate_non_empty_tissue(sim_ptr->tissue());

  const auto& species_id = sim_ptr->tissue().get_species(species_name).get_id();

  RTest<RACES::Mutants::Evolutions::SpeciesCountTest> ending_test{species_id, num_of_cells};

  RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

  sim_ptr->run(ending_test, progress_bar);
}

void SpatialSimulation::run_up_to_event(const std::string& event, const std::string& species_name,
                                        const size_t& num_of_events)
{
    run_up_to_event(event, species_name, num_of_events, false);
}

void SpatialSimulation::run_up_to_event(const std::string& event, const std::string& species_name,
                                        const size_t& num_of_events, const bool quiet)
{
  validate_non_empty_tissue(sim_ptr->tissue());

  if (event_names.count(event)==0) {
    handle_unknown_event(event);
  }

  namespace RS = RACES::Mutants::Evolutions;

  const auto& species_id = sim_ptr->tissue().get_species(species_name).get_id();

  RTest<RS::EventCountTest> ending_test{event_names.at(event), species_id, num_of_events};

  RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

  sim_ptr->run(ending_test, progress_bar);
}

void SpatialSimulation::run_until(const Logics::Formula& formula)
{
  run_until(formula, false);
}

void SpatialSimulation::run_until(const Logics::Formula& formula, const bool quiet)
{
  validate_non_empty_tissue(sim_ptr->tissue());

  RTest<RACES::Mutants::Evolutions::FormulaTest> ending_test{formula};

  RACES::UI::ProgressBar progress_bar(Rcpp::Rcout, quiet);

  sim_ptr->run(ending_test, progress_bar);
}

Rcpp::List SpatialSimulation::get_firings() const
{
  using namespace Rcpp;

  const auto last_time_sample = sim_ptr->get_statistics().get_last_time_in_history();

  auto df = get_firing_history(last_time_sample, last_time_sample);

  return DataFrame::create(_["event"]=df["event"], _["mutant"]=df["mutant"],
                           _["epistate"]=df["epistate"], _["fired"]=df["fired"]);
}

Rcpp::List SpatialSimulation::get_firing_history(const RACES::Time& minimum_time) const
{
  if (sim_ptr->get_statistics().get_history().size()==0) {
    return get_firing_history(0,0);
  }

  const auto last_time_sample = sim_ptr->get_statistics().get_last_time_in_history();

  return get_firing_history(minimum_time, last_time_sample);
}

size_t SpatialSimulation::count_history_sample_in(const RACES::Time& minimum_time,
                                           const RACES::Time& maximum_time) const
{
  size_t num_of_samples{0};
  const auto& history = sim_ptr->get_statistics().get_history();
  auto series_it = history.lower_bound(minimum_time);
  while (series_it != history.end()
         && series_it->first <= maximum_time) {
      ++num_of_samples;
      ++series_it;
  }

  return num_of_samples;
}

Rcpp::List SpatialSimulation::get_firing_history(const RACES::Time& minimum_time,
                                          const RACES::Time& maximum_time) const
{
  using namespace Rcpp;

  const size_t rows_per_sample = event_names.size()*sim_ptr->tissue().num_of_species();
  const size_t num_of_rows = count_history_sample_in(minimum_time, maximum_time)*rows_per_sample;

  CharacterVector events(num_of_rows), mutant_names(num_of_rows),
                  epi_states(num_of_rows);
  IntegerVector firings(num_of_rows);
  NumericVector times(num_of_rows);

  size_t i{0};
  const auto& history = sim_ptr->get_statistics().get_history();
  auto series_it = history.lower_bound(minimum_time);
  while (series_it != history.end() && series_it->first <= maximum_time) {
    const auto& time = series_it->first;
    const auto& t_stats = series_it->second;
    for (const auto& species: sim_ptr->tissue()) {
      for (const auto& [event_name, event_code]: event_names) {
        events[i] = event_name;
        mutant_names[i] = species.get_mutant_name();
        epi_states[i] = get_signature_string(species);

        const auto& species_it = t_stats.find(species.get_id());
        if (species_it != t_stats.end()) {
          firings[i] = count_events(species_it->second, event_code);
        } else {
          firings[i] = 0;
        }
        times[i] = time;
        ++i;
      }
    }
    ++series_it;
  }

  return DataFrame::create(_["event"]=events, _["mutant"]=mutant_names,
                           _["epistate"]=epi_states, _["fired"]=firings,
                           _["time"]=times);
}

Rcpp::List SpatialSimulation::get_count_history(const RACES::Time& minimum_time) const
{
  if (sim_ptr->get_statistics().get_history().size()==0) {
    return get_count_history(0,0);
  }

  const auto last_time_sample = sim_ptr->get_statistics().get_last_time_in_history();

  return get_count_history(minimum_time, last_time_sample);
}

Rcpp::List SpatialSimulation::get_count_history(const RACES::Time& minimum_time,
                                   const RACES::Time& maximum_time) const
{
  using namespace Rcpp;

  const size_t rows_per_sample = sim_ptr->tissue().num_of_species();
  const size_t num_of_rows = count_history_sample_in(minimum_time, maximum_time)*rows_per_sample;

  CharacterVector mutant_names(num_of_rows), epi_states(num_of_rows);
  IntegerVector counts(num_of_rows);
  NumericVector times(num_of_rows);

  size_t i{0};
  const auto& history = sim_ptr->get_statistics().get_history();
  auto series_it = history.lower_bound(minimum_time);
  while (series_it != history.end() && series_it->first <= maximum_time) {
    const auto& time = series_it->first;
    const auto& t_stats = series_it->second;
    for (const auto& species: sim_ptr->tissue()) {
      mutant_names[i] = species.get_mutant_name();
      epi_states[i] = get_signature_string(species);

      const auto& species_it = t_stats.find(species.get_id());
      if (species_it != t_stats.end()) {
        counts[i] = species_it->second.curr_cells;
      } else {
        counts[i] = 0;
      }
      times[i] = time;
      ++i;
    }
    ++series_it;
  }

  return DataFrame::create(_["mutant"]=mutant_names, _["epistate"]=epi_states,
                           _["count"]=counts, _["time"]=times);
}

Rcpp::IntegerVector SpatialSimulation::get_tissue_size() const
{
  auto size_vect = sim_ptr->tissue().size();

  return {size_vect[0], size_vect[1]};
}


RACES::Mutants::SpeciesId
get_switched_species(const RACES::Mutants::Evolutions::Tissue& tissue,
                     const RACES::Mutants::Evolutions::Species& species)
{
  auto mutant = tissue.get_mutant_species(species.get_mutant_id());

  for (const auto& mutant_species : mutant) {
    if (mutant_species.get_id() != species.get_id()) {
      return mutant_species.get_id();
    }
  }

  throw std::domain_error("The species \"" + species.get_name()
                          + "\" does not have an epigenetic status.");
}

Rcpp::List SpatialSimulation::get_rates(const std::string& species_name) const
{
  using namespace Rcpp;

  auto& species = sim_ptr->tissue().get_species(species_name);

  auto rates = List::create(_("growth") = species.get_rate(event_names.at("growth")),
                            _["death"] = species.get_rate(event_names.at("death")));

  if (species.get_methylation_signature().size()>0) {
    auto switched_id = get_switched_species(sim_ptr->tissue(), species);

    species.get_epigenetic_rate_to(switched_id);

    rates.push_back(species.get_epigenetic_rate_to(switched_id),"switch");
  }

  return rates;
}

void SpatialSimulation::update_rates(const std::string& species_name, const Rcpp::List& rates)
{
  using namespace Rcpp;
  using namespace RACES::Mutants;

  auto& species = sim_ptr->tissue().get_species(species_name);

  if (!rates.hasAttribute("names")) {
    throw std::domain_error("update_rates: The second parameter must be a Rcpp::List "
                            "with the names attribute");
  }

  auto& species_update = rate_update_history[sim_ptr->get_time()][species.get_id()];

  CharacterVector nv = rates.names();
  for (int i=0; i<nv.size(); i++) {
    auto event_name = as<std::string>(nv[i]);
    auto event_it = event_names.find(event_name);
    if (event_it == event_names.end()) {
      handle_unknown_event(event_name);
    }
    double rate = as<double>(rates[i]);
    if (event_it->second == RACES::Mutants::CellEventType::EPIGENETIC_SWITCH) {
      auto switched_id = get_switched_species(sim_ptr->tissue(), species);

      species.set_epigenetic_rate_to(switched_id, rate);
    } else {
      species.set_rate(event_it->second, rate);
    }

    species_update[event_name] = rate;
  }

  RACES::Archive::Binary::Out ruh_archive(get_rates_update_history_path());

  ruh_archive & rate_update_history;
}

Rcpp::List
SpatialSimulation::choose_cell_in(const std::string& mutant_name,
                           const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                           const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner)
{
  namespace RS = RACES::Mutants::Evolutions;

  if (sim_ptr->duplicate_internal_cells) {
    const auto rectangle = get_rectangle(lower_corner, upper_corner);
    const auto& cell = sim_ptr->choose_cell_in(mutant_name, rectangle,
                                               RACES::Mutants::CellEventType::DUPLICATION);

    return wrap_a_cell(cell);
  }

  return choose_border_cell_in(mutant_name, lower_corner, upper_corner);
}

Rcpp::List SpatialSimulation::choose_cell_in(const std::string& mutant_name)
{
  namespace RS = RACES::Mutants::Evolutions;

  if (sim_ptr->duplicate_internal_cells) {
    const auto& cell = sim_ptr->choose_cell_in(mutant_name,
                                                RACES::Mutants::CellEventType::DUPLICATION);
    return wrap_a_cell(cell);
  }

  return choose_border_cell_in(mutant_name);
}

Rcpp::List SpatialSimulation::choose_border_cell_in(const std::string& mutant_name)
{
    const auto& cell = sim_ptr->choose_border_cell_in(mutant_name);

    return wrap_a_cell(cell);
}

Rcpp::List SpatialSimulation::choose_border_cell_in(const std::string& mutant_name,
                                             const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                                             const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner)
{
    const auto rectangle = get_rectangle(lower_corner, upper_corner);
    const auto& cell = sim_ptr->choose_border_cell_in(mutant_name, rectangle);

    return wrap_a_cell(cell);
}

void SpatialSimulation::mutate_progeny(const RACES::Mutants::Evolutions::AxisPosition& x,
                                const RACES::Mutants::Evolutions::AxisPosition& y,
                                const std::string& mutated_mutant)
{
  auto pos_in_tissue = get_position_in_tissue({x,y});

  namespace RS = RACES::Mutants::Evolutions;

  sim_ptr->simulate_mutation(pos_in_tissue, mutated_mutant);
}

void SpatialSimulation::mutate_progeny(const Rcpp::List& cell_position,
                                const std::string& mutated_mutant)
{
  using namespace Rcpp;

  namespace RS = RACES::Mutants::Evolutions;

  std::vector<RACES::Mutants::Evolutions::AxisPosition> vector_position;

  for (const std::string axis : {"x", "y"}) {
    auto field = "position_"+axis;
    if (!cell_position.containsElementNamed(field.c_str())) {
      std::string msg = "Missing \"" + field + "\" element from the Rcpp::List.";

      ::Rf_error("%s", msg.c_str());
    }
    vector_position.push_back(as<RS::AxisPosition>(cell_position[field]));
  }

  return mutate_progeny(vector_position[0], vector_position[1], mutated_mutant);
}

void SpatialSimulation::sample_cells(const std::string& sample_name,
                              const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                              const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner) const
{
  using namespace RACES::Mutants;

  if (sample_name == "normal_sample") {
    ::Rf_error("Sample name \"normal_sample\" is reserved.");
  }

  auto bounding_box = get_rectangle(lower_corner, upper_corner);
  auto num_of_cells = bounding_box.size();

  RACES::Mutants::Evolutions::SampleSpecification spec(sample_name, bounding_box, num_of_cells);

  sim_ptr->sample_tissue(spec);
}

void SpatialSimulation::sample_cells(const std::string& sample_name,
                              const size_t& num_of_cells) const
{
    std::vector<RACES::Mutants::Evolutions::AxisPosition> lower_corner, upper_corner;

    for (const auto axis_size : sim_ptr->tissue().size()) {
        lower_corner.push_back(0);
        upper_corner.push_back(axis_size-1);
    }

    sample_cells(sample_name, lower_corner, upper_corner, num_of_cells);
}

void SpatialSimulation::sample_cells(const std::string& sample_name,
                              const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                              const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                              const size_t& num_of_cells) const
{
  using namespace RACES::Mutants;

  if (sample_name == "normal_sample") {
    ::Rf_error("Sample name \"normal_sample\" is reserved.");
  }

  auto bounding_box = get_rectangle(lower_corner, upper_corner);

  RACES::Mutants::Evolutions::SampleSpecification spec(sample_name, bounding_box, num_of_cells);

  sim_ptr->sample_tissue(spec);
}

SamplesForest SpatialSimulation::get_samples_forest() const
{
  return SamplesForest(*sim_ptr);
}

std::map<RACES::Mutants::SpeciesId, size_t>
count_cells_in(const RACES::Mutants::Evolutions::Tissue& tissue,
               const uint16_t& init_x, const uint16_t& init_y,
               const uint16_t& width, const uint16_t& height)
{
  std::map<RACES::Mutants::SpeciesId, size_t> counter;
  auto sizes = tissue.size();

  uint16_t x_max = std::min(static_cast<uint16_t>(init_x+width), sizes[0]);
  uint16_t y_max = std::min(static_cast<uint16_t>(init_y+height), sizes[1]);

  for (uint16_t x=init_x; x<x_max; ++x) {
    for (uint16_t y=init_y; y<y_max; ++y) {
      auto cell_proxy = tissue({x,y});
      if (!cell_proxy.is_wild_type()) {
        const RACES::Mutants::Evolutions::CellInTissue& cell = cell_proxy;

        const auto species_id = cell.get_species_id();

        auto found = counter.find(species_id);

        if (found == counter.end()) {
          counter.insert({species_id, 1});
        } else {
          ++(found->second);
        }
      }
    }
  }

  return counter;
}

inline
std::map<RACES::Mutants::SpeciesId, size_t>
count_cells_in(const RACES::Mutants::Evolutions::Tissue& tissue,
               const TissueRectangle& tumour_bounding_box,
               const uint16_t& grid_x, const uint16_t& grid_y,
               const uint16_t& width, const uint16_t& height)
{
  const uint16_t x = grid_x*width+tumour_bounding_box.lower_corner.x,
                 y = grid_y*height+tumour_bounding_box.lower_corner.y;

  return count_cells_in(tissue, x, y, width, height);
}

inline
TissueRectangle get_tissue_rectangle(const TissueRectangle& tumour_bounding_box,
                                     const uint16_t& grid_x, const uint16_t& grid_y,
                                     const uint16_t& width, const uint16_t& height)
{
  using namespace RACES::Mutants::Evolutions;

  const uint16_t x = grid_x*width+tumour_bounding_box.lower_corner.x,
                 y = grid_y*height+tumour_bounding_box.lower_corner.y;

  return TissueRectangle(PositionInTissue{x,y}, width, height);
}

std::set<RACES::Mutants::SpeciesId>
collect_species_of(const RACES::Mutants::Evolutions::Simulation& simulation,
                   const std::string& mutant_name)
{
  std::set<RACES::Mutants::SpeciesId> species_ids;

  if (mutant_name.back()!='+' && mutant_name.back()!='-') {
    const auto& tissue = simulation.tissue();

    auto mutant_id = simulation.find_mutant_id(mutant_name);

    for (const auto& species: tissue) {
      if (species.get_mutant_id()==mutant_id) {
        species_ids.insert(species.get_id());
      }
    }
  }

  return species_ids;
}

TissueRectangle SpatialSimulation::get_tumour_bounding_box() const
{
  using namespace RACES::Mutants::Evolutions;
  const auto& tissue = sim_ptr->tissue();
  const auto tissue_sizes = tissue.size();

  PositionInTissue lower_corner(static_cast<AxisSize>(tissue_sizes[0]),
                                static_cast<AxisSize>(tissue_sizes[1])), upper_corner{0,0};

  for (uint16_t grid_x=0; grid_x<tissue_sizes[0]; ++grid_x) {
    for (uint16_t grid_y=0; grid_y<tissue_sizes[1]; ++grid_y) {
      RACES::Mutants::Evolutions::PositionInTissue pos{grid_x, grid_y};
      if (!tissue(pos).is_wild_type()) {
        if (grid_x < lower_corner.x) {
          lower_corner.x = grid_x;
        }
        if (grid_y < lower_corner.y) {
          lower_corner.y = grid_y;
        }
        if (grid_x > upper_corner.x) {
          upper_corner.x = grid_x;
        }
        if (grid_y > upper_corner.y) {
          upper_corner.y = grid_y;
        }
      }
    }
  }

  return {lower_corner, upper_corner};
}

template<typename T>
inline T div_ceil(const T& dividend, const T& divisor)
{
  if (dividend==0) return dividend;

  return 1+(dividend-1)/divisor;
}

struct SpeciesConstraint
{
  const std::string species_name;
  const RACES::Mutants::SpeciesId species_id;
  const size_t min_num_of_cells;

  SpeciesConstraint(const std::string& species_name,
                   const RACES::Mutants::SpeciesId& species_id,
                   const size_t min_num_of_cells):
    species_name(species_name), species_id(species_id),
    min_num_of_cells(min_num_of_cells)
  {}

  bool is_satified(const std::map<RACES::Mutants::SpeciesId, size_t>& num_of_cells) const
  {
    auto found = num_of_cells.find(species_id);

    if (found != num_of_cells.end()) {
      return found->second >= min_num_of_cells;
    }

    return false;
  }
};

std::list<SpeciesConstraint>
get_species_constraints(const RACES::Mutants::Evolutions::Simulation& simulation,
                        const Rcpp::IntegerVector& minimum_cell_vector)
{
  std::list<SpeciesConstraint> species_constraints;

  const auto& tissue = simulation.tissue();

  const Rcpp::CharacterVector names = minimum_cell_vector.names();
  for (auto i=0; i<minimum_cell_vector.size(); ++i) {
    const auto name = Rcpp::as<std::string>(names[i]);

    if (name.back()=='+' || name.back()=='-') {
      bool name_found{false};
      for (const auto& species: tissue) {
        if (species.get_name()==name) {
          name_found = true;
          const auto& threshold = minimum_cell_vector[i];

          if (minimum_cell_vector[i] < 0) {
            throw std::domain_error("The minimum number of cells must be "
                                    "a non-negative number. Specified "
                                    + std::to_string(threshold)
                                    + " for species \"" + name + "\".");
          }
          species_constraints.emplace_back(name, species.get_id(),
                                          static_cast<size_t>(threshold));
        }
      }

      if (!name_found) {
        throw std::out_of_range("Unknown species \"" + name + "\"");
      }
    } else {
      simulation.find_mutant_id(name);
    }
  }

  return species_constraints;
}

struct MutantConstraint
{
  const std::string mutant_name;
  const std::set<RACES::Mutants::SpeciesId> species_ids;
  const size_t min_num_of_cells;

  MutantConstraint(const std::string& mutant_name,
                   const std::set<RACES::Mutants::SpeciesId>& species_ids,
                   const size_t& min_num_of_cells):
    mutant_name(mutant_name), species_ids(species_ids),
    min_num_of_cells(min_num_of_cells)
  {}

  MutantConstraint(const std::string& mutant_name,
                  std::set<RACES::Mutants::SpeciesId>&& species_ids,
                  size_t&& min_num_of_cells):
    mutant_name(mutant_name), species_ids(std::move(species_ids)),
    min_num_of_cells(min_num_of_cells)
  {}

  bool is_satified(const std::map<RACES::Mutants::SpeciesId, size_t>& num_of_cells) const
  {
    size_t total{0};

    for (const auto& species_id : species_ids) {
      auto found = num_of_cells.find(species_id);

      if (found != num_of_cells.end()) {
        total += found->second;
      }
    }

    return total > min_num_of_cells;
  }
};

std::list<MutantConstraint>
get_mutant_constraints(const RACES::Mutants::Evolutions::Simulation& simulation,
                       const Rcpp::IntegerVector& minimum_cell_vector)
{
  std::list<MutantConstraint> mutant_constraints;

  const Rcpp::CharacterVector names = minimum_cell_vector.names();
  for (auto i=0; i<minimum_cell_vector.size(); ++i) {
    const std::string name = Rcpp::as<std::string>(names[i]);

    auto species_ids = collect_species_of(simulation, name);

    if (species_ids.size()>0) {
      const auto& threshold = minimum_cell_vector[i];

      if (threshold < 0) {
        throw std::domain_error("The minimum number of cells must be "
                                "a non-negative number. Specified "
                                + std::to_string(threshold)
                                + " for species \"" + name + "\".");
      }
      mutant_constraints.emplace_back(name, std::move(species_ids),
                                      static_cast<size_t>(threshold));
    }
  }

  return mutant_constraints;
}

bool constraints_satisfied(const std::map<RACES::Mutants::SpeciesId, size_t>& cell_counts,
                           const std::list<SpeciesConstraint>& species_constraints,
                           const std::list<MutantConstraint>& mutant_constraints)
{
  for (const auto& constraint : species_constraints) {
    if (!constraint.is_satified(cell_counts)) {
      return false;
    }
  }

  for (const auto& constraint : mutant_constraints) {
    if (!constraint.is_satified(cell_counts)) {
      return false;
    }
  }

  return true;
}

std::vector<TissueRectangle>
SpatialSimulation::find_all_samples(const Rcpp::IntegerVector& minimum_cell_vector,
                             const uint16_t& width, const uint16_t& height) const
{
  auto species_constraints = get_species_constraints(*sim_ptr, minimum_cell_vector);
  auto mutant_constraints = get_mutant_constraints(*sim_ptr, minimum_cell_vector);

  auto t_bbox = get_tumour_bounding_box();

  const auto& tissue = sim_ptr->tissue();
  const auto t_width = t_bbox.upper_corner.x - t_bbox.lower_corner.x;
  const auto t_height = t_bbox.upper_corner.y - t_bbox.lower_corner.y;

  const uint16_t grid_width = t_width/width+((t_width%width>0?1:0));
  const uint16_t grid_height = t_height/height+((t_height%height>0?1:0));

  std::vector<TissueRectangle> rectangles;
  for (size_t grid_x=0; grid_x<grid_width; ++grid_x) {
    for (size_t grid_y=0; grid_y<grid_height; ++grid_y) {
      const auto cell_counts = count_cells_in(tissue, t_bbox, grid_x, grid_y,
                                              width, height);

      if (constraints_satisfied(cell_counts, species_constraints,
                                mutant_constraints)) {
        rectangles.push_back(get_tissue_rectangle(t_bbox, grid_x, grid_y,
                                                  width, height));
      }
    }
  }

  return rectangles;
}

std::vector<TissueRectangle>
SpatialSimulation::search_samples(const Rcpp::IntegerVector& minimum_cell_vector,
                           const uint16_t& width, const uint16_t& height,
                           const size_t num_of_samples, const int seed) const
{
    auto rectangles = find_all_samples(minimum_cell_vector, width, height);
    std::mt19937 gen(seed);

    if (rectangles.size()<num_of_samples) {
        std::ostringstream oss;
        if (rectangles.size()==0) {
            oss << "No box satisfies";
        } else if (rectangles.size()==1) {
            oss << "Only 1 sample satisfies";
        } else {
            oss << "Only " << rectangles.size()
                << " samples satisfy";
        }
        throw std::runtime_error(oss.str()+" the constraints.");
    }

    std::vector<TissueRectangle> output;
    while (output.size()<num_of_samples) {
        std::uniform_int_distribution<size_t> selector(0, rectangles.size()-1);

        size_t pos = selector(gen);

        output.push_back(rectangles[pos]);
        std::swap(rectangles[pos], rectangles[rectangles.size()-1]);
        rectangles.resize(rectangles.size()-1);
    }

    return output;
}

TissueRectangle SpatialSimulation::search_sample(const Rcpp::IntegerVector& minimum_cell_vector,
                                          const uint16_t& width, const uint16_t& height) const
{
  auto species_constraints = get_species_constraints(*sim_ptr, minimum_cell_vector);
  auto mutant_constraints = get_mutant_constraints(*sim_ptr, minimum_cell_vector);

  auto t_bbox = get_tumour_bounding_box();

  const auto& tissue = sim_ptr->tissue();
  const auto t_width = t_bbox.upper_corner.x - t_bbox.lower_corner.x;
  const auto t_height = t_bbox.upper_corner.y - t_bbox.lower_corner.y;

  const uint16_t grid_width = t_width/width+((t_width%width>0?1:0));
  const uint16_t grid_height = t_height/height+((t_height%height>0?1:0));

  const uint16_t diag_size = div_ceil(std::min(grid_width, grid_height),
                                      static_cast<uint16_t>(2));

  for (uint16_t diag=0; diag<diag_size; ++diag) {
    uint16_t grid_x=diag, grid_y=diag;

    for (; grid_x<grid_width-diag; ++grid_x) {
      const auto cell_counts = count_cells_in(tissue, t_bbox, grid_x, grid_y,
                                              width, height);

      if (constraints_satisfied(cell_counts, species_constraints,
                                mutant_constraints)) {
        return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
      }
    }

    for (; grid_y<grid_height-diag; ++grid_y) {
      const auto cell_counts = count_cells_in(tissue, t_bbox, grid_x, grid_y,
                                              width, height);

      if (constraints_satisfied(cell_counts, species_constraints,
                                mutant_constraints)) {
        return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
      }
    }

    for (; grid_x>diag; --grid_x) {
      const auto cell_counts = count_cells_in(tissue, t_bbox, grid_x, grid_y,
                                              width, height);

      if (constraints_satisfied(cell_counts, species_constraints,
                                mutant_constraints)) {
        return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
      }
    }

    {
      const auto cell_counts = count_cells_in(tissue, t_bbox, grid_x, grid_y,
                                              width, height);
      if (constraints_satisfied(cell_counts, species_constraints,
                                mutant_constraints)) {
        return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
      }
    }

    for (; grid_y>diag; --grid_y) {
      const auto cell_counts = count_cells_in(tissue, t_bbox, grid_x, grid_y,
                                              width, height);
      if (constraints_satisfied(cell_counts, species_constraints,
                                mutant_constraints)) {
        return get_tissue_rectangle(t_bbox, grid_x, grid_y, width, height);
      }
    }
  }
  throw std::runtime_error("No bounding box found!");
}

Logics::Variable SpatialSimulation::get_var(const std::string& name) const
{

  if ( name == "Time") {
    return Logics::Variable(sim_ptr->get_time_variable());
  }

  auto dot_pos = name.find('.');

  if (dot_pos == std::string::npos) {
    return Logics::Variable(sim_ptr->get_cardinality_variable(name));
  }

  std::string event_name(name.substr(dot_pos+1));
  std::string species_name{name.substr(0, dot_pos)};

  const auto event_id = RACES::Mutants::Evolutions::CellEvent::get_event_id(event_name);

  return Logics::Variable(sim_ptr->get_event_variable(species_name, event_id));
}

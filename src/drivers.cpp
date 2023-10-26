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

using namespace Rcpp;

struct RCloser : public Races::Drivers::Simulation::Closer
{
  bool closing() const
  {
    try {
      Rcpp::checkUserInterrupt();
    } catch (...) {
      return true;
    }

    return false;
  }
};

class Simulation : private Races::Drivers::Simulation::Simulation
{
  static bool has_names(const List& list, std::vector<std::string> aimed_names)
  {
    if (aimed_names.size() != list.size()) {
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
    if (aimed_names.size() < list.size()) {
      return false;
    }

    CharacterVector names = wrap(list.names());

    for (size_t i=0; i<names.size(); ++i) {
      if (aimed_names.count(as<std::string>(names[i]))==0) {
        return false;
      }
    }

    return true;
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
  void set_tissue_full(const std::string& name, const uint16_t& width, const uint16_t& height)
  {
    static_cast<Races::Drivers::Simulation::Simulation*>(this)->set_tissue(name, {width, height});
  }

  inline
  void set_tissue(const uint16_t& width, const uint16_t& height)
  {
    static_cast<Races::Drivers::Simulation::Simulation*>(this)->set_tissue("A tissue", {width, height});
  }

  void add_species_epigenetic(const std::string& name, const List& epigenetic_rates, 
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

    for (const std::string& states: {"+","-"}) {
      if (growth_rates.containsElementNamed(states.c_str())) {
        genotype[states].set_rate(CellEventType::DUPLICATE, as<double>(growth_rates[states]));
      }
      if (death_rates.containsElementNamed(states.c_str())) {
        genotype[states].set_rate(CellEventType::DIE, as<double>(death_rates[states]));
      }
    }

    static_cast<Races::Drivers::Simulation::Simulation*>(this)->add_species(genotype);
  }

  void add_species(const std::string& name, const double& growth_rate, const double& death_rate)
  {
    using namespace Races::Drivers;

    Genotype genotype(name, {});

    genotype[""].set_rate(CellEventType::DUPLICATE, growth_rate);
    genotype[""].set_rate(CellEventType::DIE, death_rate);

    static_cast<Races::Drivers::Simulation::Simulation*>(this)->add_species(genotype);
  }

  inline
  void add_cell(const std::string& genotype_name, const uint16_t& x, const uint16_t& y)
  {
    this->tissue().add_cell(genotype_name, {x,y});
  }

  inline
  void run_up_to(const double& time)
  {
    Races::UI::ProgressBar bar;

    RCloser closer;

    static_cast<Races::Drivers::Simulation::Simulation*>(this)->run_up_to(time, bar, closer);
  }

  List get_genotype_names() const
  {
    auto genotypes = this->tissue().get_genotypes();

    List R_genotypes(genotypes.size());

    for (size_t i=0; i<genotypes.size(); ++i) {
      R_genotypes[i] = genotypes[i].get_name();
    }

    return R_genotypes;
  }
};

RCPP_MODULE(Drivers){
  class_<Simulation>("Simulation")
  .constructor()
  .constructor<std::string>()
  .constructor<std::string, int>()
  
  .method("set_tissue", &Simulation::set_tissue_full)
  .method("set_tissue", &Simulation::set_tissue)

  .method("add_species", &Simulation::add_species_epigenetic)
  .method("add_species", &Simulation::add_species)

  .method("get_genotype_names", &Simulation::get_genotype_names)

  .method("add_cell", &Simulation::add_cell)

  .method("run_up_to", &Simulation::run_up_to);
}

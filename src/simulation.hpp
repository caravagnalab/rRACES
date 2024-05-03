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

#ifndef __RRACES_SIMULATION__
#define __RRACES_SIMULATION__

#include <vector>
#include <set>
#include <string>

#include <Rcpp.h>

#include <simulation.hpp>

#include "logics_impl.hpp"
#include "tissue_rectangle.hpp"
#include "samples_forest.hpp"

struct PlainChooser
{
  std::shared_ptr<Races::Mutants::Evolutions::Simulation> sim_ptr;
  std::string mutant_name;

  PlainChooser(const std::shared_ptr<Races::Mutants::Evolutions::Simulation>& sim_ptr,
               const std::string& mutant_name);

  inline const Races::Mutants::Evolutions::CellInTissue& operator()()
  {
    return sim_ptr->choose_cell_in(mutant_name,
                                   Races::Mutants::CellEventType::DUPLICATION);
  }
};

struct RectangularChooser : public PlainChooser
{
  Races::Mutants::RectangleSet rectangle;

  RectangularChooser(const std::shared_ptr<Races::Mutants::Evolutions::Simulation>& sim_ptr,
                     const std::string& mutant_name,
                     const std::vector<Races::Mutants::Evolutions::AxisPosition>& lower_corner,
                     const std::vector<Races::Mutants::Evolutions::AxisPosition>& upper_corner);

  inline const Races::Mutants::Evolutions::CellInTissue& operator()()
  {
    return sim_ptr->choose_cell_in(mutant_name, rectangle,
                                   Races::Mutants::CellEventType::DUPLICATION);
  }
};


class Simulation
{
  std::shared_ptr<Races::Mutants::Evolutions::Simulation> sim_ptr;  //!< The pointer to a RACES simulation object
  std::string name;      //!< The simulation name
  bool save_snapshots;   //!< A flag to preserve binary dump after object destruction

  void init(const SEXP& sexp);

  static bool has_names(const Rcpp::List& list, std::vector<std::string> aimed_names);

  static bool has_names_in(const Rcpp::List& list, std::set<std::string> aimed_names);

  static std::vector<Races::Mutants::Evolutions::Direction> get_possible_directions();

  Rcpp::List get_cells(const std::vector<Races::Mutants::Evolutions::AxisPosition>& lower_corner,
                 const std::vector<Races::Mutants::Evolutions::AxisPosition>& upper_corner,
                 const std::set<Races::Mutants::SpeciesId> &species_filter,
                 const std::set<std::string> &epigenetic_filter) const;

  Rcpp::List wrap_a_cell(const Races::Mutants::Evolutions::CellInTissue& cell) const;

  std::vector<TissueRectangle>
  find_all_samples(const Rcpp::IntegerVector& minimum_cell_vector,
                   const uint16_t& width, const uint16_t& height) const;

public:
  Simulation();

  Simulation(const SEXP& sexp);

  Simulation(const SEXP& first_param, const SEXP& second_param);

  Simulation(const std::string& simulation_name, const int& seed, const bool& save_snapshots);

  ~Simulation();

  inline void update_tissue(const std::string& name,
                            const Races::Mutants::Evolutions::AxisSize& width,
                            const Races::Mutants::Evolutions::AxisSize& height)
  {
    auto history_delta = get_history_delta();
  
    sim_ptr->set_tissue(name, {width, height});

    set_history_delta(history_delta);
  }

  inline void update_tissue(const Races::Mutants::Evolutions::AxisSize& width,
                            const Races::Mutants::Evolutions::AxisSize& height)
  {
    update_tissue(get_tissue_name(), width, height);
  }

  void add_mutant(const std::string& mutant, const Rcpp::List& epigenetic_rates,
                    const Rcpp::List& growth_rates, const Rcpp::List& death_rates);

  void add_mutant(const std::string& mutant, const double& growth_rate, const double& death_rate);

  Rcpp::List add_mutant(const Rcpp::List& prova);

  inline Races::Time get_clock() const
  {
    return sim_ptr->get_time();
  }

  void place_cell(const std::string& species_name,
                  const Races::Mutants::Evolutions::AxisPosition& x,
                  const Races::Mutants::Evolutions::AxisPosition& y);

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

  Rcpp::List get_cell(const Races::Mutants::Evolutions::AxisPosition& x,
                const Races::Mutants::Evolutions::AxisPosition& y) const;

  Rcpp::List get_cells(const std::vector<Races::Mutants::Evolutions::AxisPosition>& lower_corner,
                 const std::vector<Races::Mutants::Evolutions::AxisPosition>& upper_corner) const;

  Rcpp::List get_cells(const SEXP& first_param, const SEXP& second_param) const;

  Rcpp::List get_cells(const std::vector<std::string>& species_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  Rcpp::List get_cells(const std::vector<Races::Mutants::Evolutions::AxisPosition>& lower_corner,
                 const std::vector<Races::Mutants::Evolutions::AxisPosition>& upper_corner,
                 const std::vector<std::string>& mutant_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  Rcpp::List get_lineage_graph() const;

  Rcpp::List get_samples_info() const;

  inline void schedule_mutation(const std::string& src, const std::string& dst,
                                         const Races::Time& time)
  {
    sim_ptr->schedule_mutation(src, dst, time);
  }

  void run_up_to_time(const Races::Time& time);

  void run_up_to_size(const std::string& species_name, const size_t& num_of_cells);

  void run_up_to_event(const std::string& event, const std::string& species_name,
                       const size_t& num_of_events);

  void run_until(const Logics::Formula& formula);

  void sample_cells(const std::string& sample_name,
                    const std::vector<Races::Mutants::Evolutions::AxisPosition>& lower_corner,
                    const std::vector<Races::Mutants::Evolutions::AxisPosition>& upper_corner) const;

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
    namespace RS = Races::Mutants::Evolutions;

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

  Rcpp::List choose_border_cell_in(const std::string& mutant_name);

  Rcpp::List choose_border_cell_in(const std::string& mutant_name,
                             const std::vector<Races::Mutants::Evolutions::AxisPosition>& lower_corner,
                             const std::vector<Races::Mutants::Evolutions::AxisPosition>& upper_corner);

  Rcpp::List choose_cell_in(const std::string& mutant_name);

  Rcpp::List choose_cell_in(const std::string& mutant_name,
                      const std::vector<Races::Mutants::Evolutions::AxisPosition>& lower_corner,
                      const std::vector<Races::Mutants::Evolutions::AxisPosition>& upper_corner);

  void mutate_progeny(const Races::Mutants::Evolutions::AxisPosition& x,
                      const Races::Mutants::Evolutions::AxisPosition& y,
                      const std::string& mutated_mutant);

  void mutate_progeny(const Rcpp::List& cell_position, const std::string& mutated_mutant);

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

  TissueRectangle get_tumor_bounding_box() const;

  TissueRectangle search_sample(const Rcpp::IntegerVector& minimum_cell_vector,
                                const uint16_t& width, const uint16_t& height) const;

  inline std::vector<TissueRectangle>
  search_samples(const Rcpp::IntegerVector& minimum_cell_vector,
                 const uint16_t& width, const uint16_t& height,
                 const size_t num_of_samples) const
  {
    return search_samples(minimum_cell_vector, width, height,
                          num_of_samples, 0);
  }

  std::vector<TissueRectangle>
  search_samples(const Rcpp::IntegerVector& minimum_cell_vector,
                 const uint16_t& width, const uint16_t& height,
                 const size_t num_of_samples, int seed) const;

  Logics::Variable get_var(const std::string& name) const;
};

RCPP_EXPOSED_CLASS(Simulation)

#endif // __RRACES_SIMULATION__

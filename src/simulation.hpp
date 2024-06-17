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
  std::shared_ptr<RACES::Mutants::Evolutions::Simulation> sim_ptr;
  std::string mutant_name;

  PlainChooser(const std::shared_ptr<RACES::Mutants::Evolutions::Simulation>& sim_ptr,
               const std::string& mutant_name);

  inline const RACES::Mutants::Evolutions::CellInTissue& operator()()
  {
    return sim_ptr->choose_cell_in(mutant_name,
                                   RACES::Mutants::CellEventType::DUPLICATION);
  }
};

struct RectangularChooser : public PlainChooser
{
  RACES::Mutants::RectangleSet rectangle;

  RectangularChooser(const std::shared_ptr<RACES::Mutants::Evolutions::Simulation>& sim_ptr,
                     const std::string& mutant_name,
                     const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                     const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner);

  inline const RACES::Mutants::Evolutions::CellInTissue& operator()()
  {
    return sim_ptr->choose_cell_in(mutant_name, rectangle,
                                   RACES::Mutants::CellEventType::DUPLICATION);
  }
};


class SpatialSimulation
{
  std::shared_ptr<RACES::Mutants::Evolutions::Simulation> sim_ptr;  //!< The pointer to a RACES simulation object
  std::string name;      //!< The simulation name
  bool save_snapshots;   //!< A flag to preserve binary dump after object destruction

  using EventRateUpdateMap = std::map<std::string, double>;
  using SpeciesRateUpdateMap = std::map<RACES::Mutants::SpeciesId, EventRateUpdateMap>;
  using RateUpdateHistory = std::map<RACES::Time, SpeciesRateUpdateMap>;

  RateUpdateHistory rate_update_history;

  void init(const SEXP& sexp);

  static bool has_names(const Rcpp::List& list, std::vector<std::string> aimed_names);

  static bool has_names_in(const Rcpp::List& list, std::set<std::string> aimed_names);

  static std::vector<RACES::Mutants::Evolutions::Direction> get_possible_directions();

  Rcpp::List get_cells(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                 const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                 const std::set<RACES::Mutants::SpeciesId> &species_filter,
                 const std::set<std::string> &epigenetic_filter) const;

  Rcpp::List wrap_a_cell(const RACES::Mutants::Evolutions::CellInTissue& cell) const;

  std::vector<TissueRectangle>
  find_all_samples(const Rcpp::IntegerVector& minimum_cell_vector,
                   const uint16_t& width, const uint16_t& height) const;

  void add_mutant_rate_history(const RACES::Mutants::MutantProperties& mutant_propeties);

  inline static std::string get_rates_update_history_file_name()
  {
    return "rates_update_history.dat";
  }

  inline std::filesystem::path get_rates_update_history_path() const
  {
    return std::filesystem::absolute(sim_ptr->get_logger().get_directory()/
                                        get_rates_update_history_file_name());
  }
public:

  template<typename SAMPLES>
  static Rcpp::List get_samples_info(const SAMPLES& samples)
  {
    using namespace Rcpp;
    CharacterVector sample_name(samples.size());
    NumericVector time(samples.size());
    IntegerVector ymin(samples.size()), ymax(samples.size()),
                  xmin(samples.size()), xmax(samples.size()),
                  num_of_cells(samples.size()),
                  tumour_in_bbox(samples.size());

    size_t i{0};
    for (const auto& sample : samples) {
        sample_name[i] = sample.get_name();
        time[i] = sample.get_time();
        num_of_cells[i] = sample.get_cell_ids().size();
        tumour_in_bbox[i] = sample.get_tumour_cells_in_bbox();
        const auto& bounding_box = sample.get_bounding_box();
        xmin[i] = bounding_box.lower_corner.x;
        xmax[i] = bounding_box.upper_corner.x;
        ymin[i] = bounding_box.lower_corner.y;
        ymax[i] = bounding_box.upper_corner.y;
        ++i;
    }

    return DataFrame::create(_["name"]=sample_name, _["xmin"]=xmin,
                             _["ymin"]=ymin, _["xmax"]=xmax,
                             _["ymax"]=ymax,
                             _["tumour_cells"]=num_of_cells,
                             _["tumour_cells_in_bbox"]=tumour_in_bbox,
                             _["time"]=time);
  }

  SpatialSimulation();

  SpatialSimulation(const SEXP& sexp);

  SpatialSimulation(const SEXP& first_param, const SEXP& second_param);

  SpatialSimulation(const std::string& simulation_name, const SEXP& seed,
                    const bool& save_snapshots);

  SpatialSimulation(const std::string& simulation_name, const bool& seed,
                    const bool& save_snapshots);

  ~SpatialSimulation();

  inline void update_tissue(const std::string& name,
                            const RACES::Mutants::Evolutions::AxisSize& width,
                            const RACES::Mutants::Evolutions::AxisSize& height)
  {
    auto history_delta = get_history_delta();
  
    sim_ptr->set_tissue(name, {width, height});

    set_history_delta(history_delta);
  }

  inline void update_tissue(const RACES::Mutants::Evolutions::AxisSize& width,
                            const RACES::Mutants::Evolutions::AxisSize& height)
  {
    update_tissue(get_tissue_name(), width, height);
  }

  void add_mutant(const std::string& mutant, const Rcpp::List& epigenetic_rates,
                  const Rcpp::List& growth_rates, const Rcpp::List& death_rates);

  void add_mutant(const std::string& mutant, const double& growth_rate,
                  const double& death_rate);

  inline RACES::Time get_clock() const
  {
    return sim_ptr->get_time();
  }

  void place_cell(const std::string& species_name,
                  const RACES::Mutants::Evolutions::AxisPosition& x,
                  const RACES::Mutants::Evolutions::AxisPosition& y);

  size_t count_history_sample_in(const RACES::Time& minimum_time,
                                 const RACES::Time& maximum_time) const;

  Rcpp::List get_added_cells() const;

  Rcpp::List get_counts() const;

  inline Rcpp::List get_count_history() const
  {
    return get_count_history(0);
  }

  Rcpp::List get_count_history(const RACES::Time& minimum_time) const;

  Rcpp::List get_count_history(const RACES::Time& minimum_time,
                         const RACES::Time& maximum_time) const;

  Rcpp::List get_cells() const;

  Rcpp::List get_cell(const RACES::Mutants::Evolutions::AxisPosition& x,
                const RACES::Mutants::Evolutions::AxisPosition& y) const;

  Rcpp::List get_cells(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                 const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner) const;

  Rcpp::List get_cells(const SEXP& first_param, const SEXP& second_param) const;

  Rcpp::List get_cells(const std::vector<std::string>& species_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  Rcpp::List get_cells(const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                 const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                 const std::vector<std::string>& mutant_filter,
                 const std::vector<std::string>& epigenetic_filter) const;

  Rcpp::List get_lineage_graph() const;

  inline Rcpp::List get_samples_info() const
  {
    return get_samples_info(sim_ptr->get_tissue_samples());
  }

  inline void schedule_mutation(const std::string& src, const std::string& dst,
                                         const RACES::Time& time)
  {
    sim_ptr->schedule_mutation(src, dst, time);
  }

  void run_up_to_time(const RACES::Time& time);

  void run_up_to_size(const std::string& species_name, const size_t& num_of_cells);

  void run_up_to_event(const std::string& event, const std::string& species_name,
                       const size_t& num_of_events);

  void run_until(const Logics::Formula& formula);

  void sample_cells(const std::string& sample_name, const size_t& num_of_cells) const;

  void sample_cells(const std::string& sample_name,
                    const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                    const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner) const;

  void sample_cells(const std::string& sample_name,
                    const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                    const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner,
                    const size_t& num_of_cells) const;

  Rcpp::List get_firings() const;

  inline Rcpp::List get_firing_history() const
  {
    return get_firing_history(0);
  }

  Rcpp::List get_firing_history(const RACES::Time& minimum_time) const;

  Rcpp::List get_firing_history(const RACES::Time& minimum_time,
                          const RACES::Time& maximum_time) const;

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

  Rcpp::List get_rates_update_history() const;

  void update_rates(const std::string& species_name, const Rcpp::List& list);

  template<typename CHOOSER, std::enable_if_t<std::is_base_of_v<PlainChooser, CHOOSER>, bool> = true>
  Rcpp::List choose_border_cell_in(CHOOSER& chooser)
  {
    namespace RS = RACES::Mutants::Evolutions;

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
                             const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                             const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner);

  Rcpp::List choose_cell_in(const std::string& mutant_name);

  Rcpp::List choose_cell_in(const std::string& mutant_name,
                      const std::vector<RACES::Mutants::Evolutions::AxisPosition>& lower_corner,
                      const std::vector<RACES::Mutants::Evolutions::AxisPosition>& upper_corner);

  void mutate_progeny(const RACES::Mutants::Evolutions::AxisPosition& x,
                      const RACES::Mutants::Evolutions::AxisPosition& y,
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

  inline bool is_border_growth_model() const
  {
    return !sim_ptr->duplicate_internal_cells;
  }

  inline void set_border_growth_model(const bool border_growth_model)
  {
    sim_ptr->duplicate_internal_cells = !border_growth_model;
  }

  inline RACES::Time get_history_delta() const
  {
    return sim_ptr->get_statistics().get_history_delta();
  }

  inline void set_history_delta(const RACES::Time history_time_delta)
  {
    sim_ptr->get_statistics().set_history_delta(history_time_delta);
  }

  static SpatialSimulation load(const std::string& directory_name);

  SamplesForest get_samples_forest() const;

  TissueRectangle get_tumour_bounding_box() const;

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
  
  static SpatialSimulation build_simulation(const SEXP& simulation_name, const SEXP& width,
                                            const SEXP& height, const SEXP& save_snapshots,
                                            const SEXP& seed);
};

RCPP_EXPOSED_CLASS(SpatialSimulation)

#endif // __RRACES_SIMULATION__

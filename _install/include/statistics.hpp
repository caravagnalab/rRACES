/**
 * @file statistics.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines simulation statistics
 * @version 0.9
 * @date 2023-10-12
 * 
 * @copyright Copyright (c) 2023
 * 
 * MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __RACES_STATISTICS__
#define __RACES_STATISTICS__

#include <list>

#include "species.hpp"
#include "cell_event.hpp"


namespace Races 
{

namespace Drivers
{

namespace Simulation
{

class TissueStatistics;

/**
 * @brief A class for accounting simulation species statistics
 */
struct SpeciesStatistics 
{

    Time rise_time;         //!< the time of the first appearance 
    Time extinction_time;   //!< the time of the last appearance

    size_t total_cells;     //!< the number of cells that have been in the species
    size_t curr_cells;      //!< the number of cells currently in the species
    size_t killed_cells;    //!< the number of species cells that were killed
    size_t lost_cells;      //!< the number of cells that have overcome the tissue border

    /**
     * @brief The empty constructor
     */
    SpeciesStatistics();

    /**
     * @brief A constructor
     * 
     * @param num_of_cells is the number of cells in the species
     */
    explicit SpeciesStatistics(const size_t& num_of_cells);

    /**
     * @brief Save species statistics in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & rise_time
                & extinction_time
                & total_cells
                & curr_cells
                & killed_cells
                & lost_cells;
    }
    
    /**
     * @brief Load species statistics from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded species statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SpeciesStatistics load(ARCHIVE& archive)
    {
        SpeciesStatistics stats;

        archive & stats.rise_time
                & stats.extinction_time
                & stats.total_cells
                & stats.curr_cells
                & stats.killed_cells
                & stats.lost_cells;

        return stats;
    }

    friend class TissueStatistics;
};

/**
 * @brief A class to collect simulation statistics about a tissue
 */
class TissueStatistics
{
    using time_point = std::chrono::steady_clock::time_point;

    std::map<EpigeneticGenotypeId, SpeciesStatistics> s_statistics;  //!< a map from species id to statistics
    
    std::list<Time> sim_times;             //!< the simulated times of the last recorded events
    std::list<time_point> real_times;      //!< the recording times of the last events 

    size_t max_stored_times;               //!< the maximum number of times to store

    size_t total_events;                   //!< number of recorded events
    time_point first_event_time;           //!< first recoded event time
public:
    /**
     * @brief An empty constructor
     */
    TissueStatistics();

    /**
     * @brief Get the statistics of a species
     * 
     * @param species is the species whose statistics are aimed
     * @return a non-constant reference to the statistics of `species` 
     */
    inline SpeciesStatistics& operator[](const Species& species)
    {
        return s_statistics.at(species.get_id());
    }

    /**
     * @brief Get the statistics of a species
     * 
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return a non-constant reference to the statistics of species having 
     *      `species_id` as identifier 
     */
    inline SpeciesStatistics& operator[](const EpigeneticGenotypeId& species_id)
    {
        return s_statistics.at(species_id);
    }

    /**
     * @brief Get the statistics of a species
     * 
     * @param species is the species whose statistics are aimed
     * @return a non-constant reference to the statistics of `species` 
     */
    inline const SpeciesStatistics& at(const Species& species) const
    {
        return s_statistics.at(species.get_id());
    }

    /**
     * @brief Get the statistics of a species
     * 
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return a constant reference to the statistics of species having 
     *      `species_id` as identifier 
     */
    inline const SpeciesStatistics& at(const EpigeneticGenotypeId& species_id) const
    {
        return s_statistics.at(species_id);
    }

    /**
     * @brief Test whether the object contains statistics for a species
     * 
     * @param species_id is the identifier of the species whose statistics are aimed
     * @return `true` if and only if the object contains statistics for the 
     *          specified species
     */
    inline bool contains_data_for(const EpigeneticGenotypeId& species_id) const
    {
        return s_statistics.count(species_id)==1;
    }

    /**
     * @brief Test whether the object contains statistics for a species
     * 
     * @param species is the species whose statistics are aimed
     * @return `true` if and only if the object contains statistics for the 
     *          specified species
     */
    bool contains_data_for(const Species& species) const; 

    /**
     * @brief Record a death
     * 
     * @param genotype_id is the genotype id of the dying cell
     * @param time is the death time
     */
    void record_death(const EpigeneticGenotypeId& genotype_id, const Time &time);

    /**
     * @brief Record a lost cell
     * 
     * A cell is lost whenever it is pushed outside the tissue border. 
     * This method record in the statistics that a cell is lost.
     * 
     * @param genotype_id is the genotype id of the lost cell
     * @param time is the time in which the lost cell has been pushed 
     *       outside the tissue border
     */
    void record_lost(const EpigeneticGenotypeId& genotype_id, const Time &time);

    /**
     * @brief Record a cell duplication
     * 
     * @param genotype_id is the driver genotype id of the duplicating cell
     */
    void record_duplication(const EpigeneticGenotypeId& genotype_id);

    /**
     * @brief Record a driver mutation
     * 
     * @param initial_id is the original driver genotype id of the mutating cell 
     * @param final_id is the driver genotype id of the mutated cell
     * @param time is the epigenetic event time
     */
    void record_mutation(const EpigeneticGenotypeId& initial_id, const EpigeneticGenotypeId& final_id, const Time &time);

    /**
     * @brief Record a cell duplication and an epigenetic event
     * 
     * @param genotype_id is the driver genotype id of the duplicating cell 
     * @param epigenetic_genotype is the driver genotype id of the mutated cell
     * @param time is the epigenetic event time
     */
    void record_duplication_epigenetic_event(const EpigeneticGenotypeId& genotype_id, const EpigeneticGenotypeId& epigenetic_genotype, const Time &time);

    /**
     * @brief Record the last event
     * 
     * This method records a new event and, whenever the 
     * queues of the recorded times overlaps the maximum size, 
     * deletes the oldest recorded times.
     * 
     * @param event is the event to record
     * @param time is the event time
     */
    void record_event(const CellEvent& event, const Time &time);

    /**
     * @brief Get the number of recorded times
     * 
     * @return the number of recorded times
     */
    const size_t& get_recorded_time_number() const;

    /**
     * @brief Get the elapsed time
     * 
     * @return the elapsed time
     */
    inline std::chrono::steady_clock::duration get_elapsed_time() const
    {
        return std::chrono::steady_clock::now()-first_event_time;
    }

    /**
     * @brief Get the simulated time
     * 
     * @return the simulated time
     */
    inline Time get_simulated_time() const
    {
        return sim_times.back();
    }

    /**
     * @brief Get the number of the last recorded events over time
     * 
     * @return the number of the last recorded events over time
     */
    template<class ToDuration> 
    double get_last_recorded_events_over_time() const
    {
        using namespace std::chrono;

        auto time = duration_cast<ToDuration>(real_times.back()-real_times.front()).count();

        return static_cast<double>(real_times.size())/time;
    }
    
    /**
     * @brief Save tissue statistics in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & s_statistics
                & sim_times
                & real_times
                & max_stored_times
                & total_events
                & first_event_time;
    }
    
    /**
     * @brief Load tissue statistics from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded tissue statistics
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static TissueStatistics load(ARCHIVE& archive)
    {
        TissueStatistics stats;

        archive & stats.s_statistics
                & stats.sim_times
                & stats.real_times
                & stats.max_stored_times
                & stats.total_events
                & stats.first_event_time;

        return stats;
    }
};

}   // Simulation

}   // Drivers

}   // Races


#endif // __RACES_STATISTICS__
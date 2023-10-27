/**
 * @file rate_update.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines liveness rate updates
 * @version 0.1
 * @date 2023-10-18
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

#ifndef __RACES_RATE_UPDATE__
#define __RACES_RATE_UPDATE__

#include "driver_genotype.hpp"
#include "cell_event.hpp"
#include "simulation_event.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief A structure to represent liveness rate update
 */
struct RateUpdate : public SimulationEvent
{
    using Type = SimulationEvent::Type;

    EpigeneticGenotypeId species_id;   //!< The involved species id
    CellEventType event_type;          //!< The liveness event type
    double new_rate;                   //!< The new rate

    /**
     * @brief A constructor
     * 
     * @param species_id is the identifier of the species whose rate is changed by the event
     * @param event_type is the event type whose rate is changed by the event
     * @param new_rate is the new rate for `event_type` in the species `species_id`
     */
    RateUpdate(const EpigeneticGenotypeId& species_id, 
               const CellEventType& event_type, const double& new_rate);

    /**
     * @brief A constructor
     * 
     * @param species is the species whose rate is changed by the event
     * @param event_type is the event type whose rate is changed by the event
     * @param new_rate is the new rate for `event_type` in the species `species_id`
     */
    RateUpdate(const EpigeneticGenotype& species, const CellEventType& event_type,
               const double& new_rate);

    /**
     * @brief Save a timed genomic mutation in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & species_id
                & event_type
                & new_rate;
    }

    /**
     * @brief Load a timed genomic mutation from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static RateUpdate load(ARCHIVE& archive)
    {
        EpigeneticGenotypeId species_id;
        CellEventType event_type;
        double new_rate;

        archive & species_id
                & event_type
                & new_rate;

        return {species_id, event_type, new_rate};
    }

    inline Type type() const {
        return Type::LIVENESS_RATE_UPDATE;
    } 
};

}   // Simulation

}   // Drivers

}   // Races


/**
 * @brief Test the equivalence between two liveness rate updates
 * 
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two liveness rate updates represent
 *      the same event
 */
inline bool operator==(const Races::Drivers::Simulation::RateUpdate& lhs, 
                       const Races::Drivers::Simulation::RateUpdate& rhs)
{
    return (lhs.species_id == rhs.species_id)
            && (lhs.event_type == rhs.event_type)
            && (lhs.new_rate == rhs.new_rate);
}

#endif // __RACES_RATE_UPDATE__
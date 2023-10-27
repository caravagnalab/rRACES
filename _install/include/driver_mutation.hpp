/**
 * @file driver_mutation.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genomic mutations
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

#ifndef __RACES_DRIVER_MUTATION__
#define __RACES_DRIVER_MUTATION__

#include "driver_genotype.hpp"
#include "simulation_event.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief A structure to represent driver genomic mutation
 */
struct DriverMutation : public SimulationEvent
{
    using Type = SimulationEvent::Type;

    GenotypeId initial_id;   //!< The initial genomic genotype identifier
    GenotypeId final_id;     //!< The final genomic genotype identifier

    /**
     * @brief A constructor
     * 
     * @param initial_id is the identifier of the initial genomic genotype
     * @param final_id is the identifier of the final genomic genotype
     */
    DriverMutation(const GenotypeId& initial_id, const GenotypeId& final_id);

    /**
     * @brief A constructor
     * 
     * @param initial_genotype is the initial genomic genotype
     * @param final_genotype is the final genomic genotype
     */
    DriverMutation(const Genotype& initial_genotype, const Genotype& final_genotype);

    /**
     * @brief Save a genomic mutation in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & initial_id 
                & final_id;
    }

    /**
     * @brief Load a genomic mutation from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static DriverMutation load(ARCHIVE& archive)
    {
        GenotypeId initial_id;
        GenotypeId final_id;

        archive & initial_id 
                & final_id;

        return {initial_id, final_id};
    }

    inline Type type() const {
        return Type::DRIVER_MUTATION;
    } 
};

}   // Simulation

}   // Drivers

}   // Races

/**
 * @brief Test the equivalence between two driver mutations
 * 
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two driver mutations represent
 *      the same event
 */
inline bool operator==(const Races::Drivers::Simulation::DriverMutation& lhs,
                       const Races::Drivers::Simulation::DriverMutation& rhs)
{
    return (lhs.initial_id == rhs.initial_id) 
            && (lhs.final_id == rhs.final_id);
}

#endif // __RACES_DRIVER_MUTATION__
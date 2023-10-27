/**
 * @file timed_event.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines timed events
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

#ifndef __RACES_TIMED_EVENT__
#define __RACES_TIMED_EVENT__

#include <map>
#include <string>
#include <functional>

#include "time.hpp"
#include "archive.hpp"
#include "event_wrapper.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief A structure to represent timed driver genomic mutation
 */
struct TimedEvent : public SimulationEventWrapper
{
    using Type = SimulationEvent::Type;

    Time time;  //!< The event time
    
    /**
     * @brief A constructor
     * 
     * @param time is the simulation time of the event
     * @param event is the event
     */
    TimedEvent(const Time& time, const SimulationEventWrapper& event);

    /**
     * @brief Save a timed event in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & time
                & static_cast<const SimulationEventWrapper&>(*this);
    }

    /**
     * @brief Load a timed event from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed event
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static TimedEvent load(ARCHIVE& archive)
    {
        Time time;

        archive & time;

        return {time, SimulationEventWrapper::load(archive)};
    }
};

}   // Simulation

}   // Drivers

}   // Races


template<>
struct std::greater<Races::Drivers::Simulation::TimedEvent> {
    inline constexpr bool operator()(const Races::Drivers::Simulation::TimedEvent &lhs, 
                                     const Races::Drivers::Simulation::TimedEvent &rhs) const 
    {
        return lhs.time > rhs.time;
    }
};

/**
 * @brief Test the equivalence between two timed events
 * 
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two timed events represent
 *      the same event
 */
inline bool operator==(const Races::Drivers::Simulation::TimedEvent& lhs, 
                       const Races::Drivers::Simulation::TimedEvent& rhs)
{
    using namespace Races::Drivers::Simulation;

    return (lhs.time == rhs.time) 
            && (static_cast<const SimulationEventWrapper&>(lhs)
                ==static_cast<const SimulationEventWrapper&>(rhs));
}

#endif // __RACES_TIMED_EVENT__
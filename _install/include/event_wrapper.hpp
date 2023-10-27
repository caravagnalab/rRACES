/**
 * @file event_wrapper.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a simulation event wrapper
 * @version 0.2
 * @date 2023-10-25
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

#ifndef __RACES_EVENT_WRAPPER__
#define __RACES_EVENT_WRAPPER__

#include <memory>

#include "archive.hpp"
#include "driver_mutation.hpp"
#include "rate_update.hpp"
#include "sampling.hpp"

namespace Races 
{

namespace Drivers
{

namespace Simulation
{

/**
 * @brief A structure to represent a generic simulation event
 */
struct SimulationEventWrapper
{
    using Type = SimulationEvent::Type;

    Type type;  //!< The type of the simulation event
    std::shared_ptr<void> event;    //!< The pointer to the real event

    /**
     * @brief A constructor that wraps a driver mutation event
     * 
     * @param driver_mutation is the event to wrap
     */
    SimulationEventWrapper(const DriverMutation& driver_mutation);

    /**
     * @brief A constructor that wraps a liveness rate event
     * 
     * @param rate_update is the event to wrap
     */
    SimulationEventWrapper(const RateUpdate& rate_update);

    /**
     * @brief A constructor that wraps a sampling event
     * 
     * @param sampling is the event to wrap
     */
    SimulationEventWrapper(const Sampling& sampling);

    /**
     * @brief Get the wrapped event
     * 
     * @tparam EVENT_TYPE is the type of the wrapped event
     * @return a reference to the wrapped event
     */
    template<typename EVENT_TYPE>
    EVENT_TYPE& get_event()
    {
        return *(std::static_pointer_cast<EVENT_TYPE>(event));
    }

    /**
     * @brief Get the wrapped event
     * 
     * @tparam EVENT_TYPE is the type of the wrapped event
     * @return a constant reference to the wrapped event 
     */
    template<typename EVENT_TYPE>
    const EVENT_TYPE& get_event() const
    {
        return *(std::static_pointer_cast<EVENT_TYPE>(event));
    }

    /**
     * @brief Save a simulation event in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & type;

        switch(type) {
            case Type::DRIVER_MUTATION:
                archive & get_event<DriverMutation>();
                break;
            case Type::LIVENESS_RATE_UPDATE:
                archive & get_event<RateUpdate>();
                break;
            case Type::SAMPLING:
                archive & get_event<Sampling>();
                break;
            default:
                throw std::runtime_error("Unsupported event type");
        }
    }

    /**
     * @brief Load a timed event from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed event
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SimulationEventWrapper load(ARCHIVE& archive)
    {
        Type type;

        archive & type;

        switch(type) {
            case Type::DRIVER_MUTATION:
                {
                    auto mutation = DriverMutation::load(archive);

                    return SimulationEventWrapper(mutation);
                }
            case Type::LIVENESS_RATE_UPDATE:
                {
                    auto rate_update = RateUpdate::load(archive);

                    return SimulationEventWrapper(rate_update);
                }
            case Type::SAMPLING:
                {
                    auto sampling = Sampling::load(archive);

                    return SimulationEventWrapper(sampling);
                }
            default:
                throw std::runtime_error("Unsupported event type");
        }
    }
};

}   // Simulation

}   // Drivers

}   // Races

/**
 * @brief Test the equivalence between two event wrappers
 * 
 * @param lhs is the left-hand side of the equivalence
 * @param rhs is the right-hand side of the equivalence
 * @return `true` if and only if the two wrappers represent
 *      the same event
 */
bool operator==(const Races::Drivers::Simulation::SimulationEventWrapper& lhs, 
                const Races::Drivers::Simulation::SimulationEventWrapper& rhs);

#endif // __RACES_EVENT_WRAPPER__
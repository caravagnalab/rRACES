/**
 * @file logger.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines simulation loggers
 * @version 0.11
 * @date 2023-10-23
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

#ifndef __RACES_LOGGER__
#define __RACES_LOGGER__

#include <string>
#include <filesystem>

#include "cell.hpp"
#include "cell_event.hpp"
#include "time.hpp"
#include "tissue_sample.hpp"

namespace Races 
{

namespace Drivers 
{

namespace Simulation 
{

class Simulation;

/**
 * @brief The simulator logger concept
 */
class BasicLogger
{
protected:
    std::filesystem::path directory;   //!< the log directory

    /**
     * @brief Get a string representing the current time
     * 
     * This method produces a string in the format
     *    {year}{month}{day}_{hour}{minutes}{seconds}
     * 
     * @return a string representing the current time 
     */
    static std::string get_time_string();
public:
    /**
     * @brief The empty constructor
     */
    BasicLogger();

    /**
     * @brief A constructor
     * 
     * @param simulation_dir is the new simulation logging directory
     */
    explicit BasicLogger(const std::filesystem::path simulation_dir);

    /**
     * @brief Record an event
     * 
     * @param type is the event type
     * @param cell is the cell on which event has been occurred
     * @param time it the event time
     */
    void record(const CellEventType& type, const CellInTissue& cell, const Time& time);

    /**
     * @brief Record an initial cell
     * 
     * @param cell is the initial cell to record 
     */
    void record_initial_cell(const CellInTissue& cell);

    /**
     * @brief Save a simulation snapshot
     * 
     * @param simulation is the simulation whose snapshot is requested
     */
    void snapshot(const Simulation& simulation);

    /**
     * @brief Flush archive data
     */
    inline void flush_archives()
    {}

    /**
     * @brief Save tissue sample 
     * 
     * @param simulation_dir is the path of the simulation directory
     * @param tissue_sample is the tissue sample to log
     */
    static void 
    save_sample(const std::filesystem::path simulation_dir,
                const Races::Drivers::Simulation::TissueSample& tissue_sample);

    /**
     * @brief Close open archives
     */
    virtual void close()
    {}

    /**
     * @brief Reset the logger
     * 
     * @param new_directory is the new logger directory
     */
    inline void reset(const std::filesystem::path new_directory)
    {
        close();

        directory = new_directory;
    }

    /**
     * @brief Rename the logger directory
     * 
     * @param new_path is the new path of the logger directory
     */
    inline void rename_directory(const std::filesystem::path& new_path)
    {
        close();

        if (std::filesystem::exists(directory)) {
            std::filesystem::rename(directory, new_path);
        }

        directory = new_path;
    }

    /**
     * @brief Get the log directory
     * 
     * @return the directory containing all the logs
     */
    inline const std::filesystem::path& get_directory() const
    {
        return directory;
    }

    /**
     * @brief Save tissue sample 
     * 
     * @param tissue_sample is the tissue sample to log
     */
    inline void save_sample(const Races::Drivers::Simulation::TissueSample& tissue_sample) const
    {
        BasicLogger::save_sample(directory, tissue_sample);
    }
};

}   // Simulation

}   // Drivers

}   // Races

#endif // __RACES_LOGGER__
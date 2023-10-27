/**
 * @file progress_bar.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a progress bar
 * @version 0.11
 * @date 2023-10-26
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

#ifndef __RACES_PROGRESS_BAR__
#define __RACES_PROGRESS_BAR__

#include <chrono>
#include <string>

#include "variables.hpp"

#if WITH_INDICATORS
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

#include "indicators.hpp"

#pragma GCC diagnostic pop
#endif // WITH_INDICATORS

namespace Races
{

namespace UI
{

/**
 * @brief A class to represent progress bars
 * 
 * This class avoids too many effective updates by 
 * setting a minimum time interval between two of 
 * them. Any update occurring less than the specified
 * time interval from the last update is recorded in 
 * the status, but not graphically visualized.
 */
class ProgressBar
{
    std::chrono::system_clock::time_point last_update;  //!< the last update time

    unsigned int percentage;  //!< the percentage level
    std::string message;      //!< the progress bar message

    bool updated;             //!< the last percentage has been shown

#if WITH_INDICATORS
    indicators::ProgressBar* indicator;  //!< the progress bar implementation
#endif

public:
    std::chrono::system_clock::duration update_interval;  //!< the time between two effective updates

    /**
     * @brief The constructor
     */
    ProgressBar();

    /**
     * @brief A constructor
     * 
     * @param hide is a hide/show flag 
     */
    ProgressBar(const bool hide);

    /**
     * @brief Set a message in the progress bar
     * 
     * @param message is the message
     * @return a reference to the updated object
     */
    ProgressBar& set_message(const std::string& message);

    /**
     * @brief Update the progress bar level
     * 
     * @param percentage is the percentage level
     * @return a reference to the updated object
     */
    ProgressBar& set_progress(const unsigned int percentage);

    /**
     * @brief Update the progress bar level
     * 
     * @param percentage is the percentage level
     * @param message is the message to be set
     * @return a reference to the updated object
     */
    ProgressBar& set_progress(const unsigned int percentage, const std::string& message);

    /**
     * @brief Get the percentage level
     * 
     * @return a constant reference to the percentage level
     */
    const unsigned int& get_progress() const;

    /**
     * @brief Update elapse time
     * 
     * @return a reference to the updated object
     */
    inline ProgressBar& update_elapsed_time()
    {
        return set_progress(get_progress());
    }

    /**
     * @brief Show the console cursor
     */
    static void show_console_cursor();

    /**
     * @brief Hide the console cursor
     */
    static void hide_console_cursor();

    /**
     * @brief The constructor
     */
    ~ProgressBar();
};

}  // UI

} // Races

#endif // __RACES_PROGRESS_BAR__
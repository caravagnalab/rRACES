/**
 * @file binary_logger.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a binary simulation logger
 * @version 0.17
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

#ifndef __RACES_BINARY_LOGGER__
#define __RACES_BINARY_LOGGER__

#include <fstream>
#include <filesystem>
#include <list>

#include "logger.hpp"

namespace Races 
{

namespace Drivers 
{

namespace Simulation 
{

/**
 * @brief A binary logger
 * 
 * The objects of this class record cell duplications and 
 * tissue status in binary files 
 */
struct BinaryLogger : public BasicLogger
{
    Archive::Binary::Out cell_archive; //!< the current cell output file
    size_t cells_per_file;             //!< the number of cells to be saved before changing file
    size_t cell_in_current_file;       //!< cells saved in current file
    uint16_t current_file_number;      //!< the current file number

    /**
     * @brief Get a new snapshot file full path
     * 
     * @return a new snapshot file full path
     */
    std::filesystem::path get_snapshot_path() const;

    /**
     * @brief Open a new cell file
     * 
     * This method close the last opened cell file and open a new one.
     */
    void rotate_cell_file();

    /**
     * @brief Record a cell
     * 
     * @param cell is the cell on which event has been occurred
     */
    void record_cell(const CellInTissue& cell);
public:
    /**
     * @brief A structure for cell recovering
     */
    class CellReader {
        std::filesystem::path directory;    //!< cell file directory
        size_t cells_per_archive;           //!< number of cells per archive
        size_t bytes_per_cell;              //!< number of bytes per cell record

        uint16_t last_file_number;          //!< last cell file number in the directory
        uint64_t number_of_cells;           //!< total number of cells in the directory

        uint16_t cell_archive_id;           //!< number of the currently open cell archive
        Archive::Binary::In cell_archive;   //!< stream of the currently open cell archive

        /**
         * @brief Get the position in the archives
         * 
         * @param cell_id 
         * @return std::streampos 
         */
        inline std::streampos get_cell_pos(const CellId& cell_id) const
        {
            return (cell_id % cells_per_archive)*bytes_per_cell;
        }

    public:
        /**
         * @brief The cell reader constructor
         * 
         * @param directory is the directory containing the cell files
         */
        explicit CellReader(const std::filesystem::path& directory);

        /**
         * @brief Get a timed-labelled cell from the directory
         * 
         * @param cell_id is the identifier of the aimed cell
         * @return the cell in the `cell_id`-th position of the 
         *      cell archives
         */
        Cell operator[](const CellId& cell_id);
    };

    /**
     * @brief Get the next cell file full path
     * 
     * @return the next cell file full path
     */
    static std::filesystem::path get_cell_archive_path(const std::filesystem::path& directory, const uint16_t& file_number);

    /**
     * @brief The empty constructor
     */
    BinaryLogger();

    /**
     * @brief A constructor
     * 
     * @param simulation_dir is the simulation logger directory name
     * @param cells_per_file is the number of cells per file
     */
    explicit BinaryLogger(const std::filesystem::path simulation_dir, const size_t cells_per_file=1<<27);

    /**
     * @brief The copy-constructor
     * 
     * @param orig is the template for the new object
     */
    BinaryLogger(const BinaryLogger& orig);

    /**
     * @brief The copy operator
     * 
     * @param orig is the template for the new object
     * @return a reference to the updated object
     */
    BinaryLogger& operator=(const BinaryLogger& orig);

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
    {
        if (cell_archive.is_open()) {
            cell_archive.flush();
        }
    }

    /**
     * @brief Close open archives
     */
    inline void close() override
    {
        if (cell_archive.is_open()) {
            cell_archive.close();
        }
    }

    /**
     * @brief Reset the logger
     * 
     * @param output_directory is the output directory name
     */
    void reset(const std::string& directory);

    /**
     * @brief Reset the logger
     */
    inline void reset()
    {
        reset(directory);
    }

    /**
     * @brief Save tissue sample 
     * 
     * @param tissue_sample is the tissue sample to log
     */
    inline void save_sample(const Races::Drivers::Simulation::TissueSample& tissue_sample) const
    {
        BinaryLogger::save_sample(directory, tissue_sample);
    }

    /**
     * @brief Load the saved samples
     * 
     * @return the list of the saved samples
     */
    inline std::list<Races::Drivers::Simulation::TissueSample>
    load_samples() const
    {
        return BinaryLogger::load_samples(directory);
    }

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
     * @brief Load the saved samples
     * 
     * @param simulation_dir is the path of the simulation directory
     * @return the list of the saved samples
     */
    static std::list<Races::Drivers::Simulation::TissueSample>
    load_samples(const std::filesystem::path simulation_dir);

    /**
     * @brief Save a `BinaryLogger` in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<std::string>(directory)
                & cells_per_file
                & cell_in_current_file
                & current_file_number;
    }

    /**
     * @brief Load a `BinaryLogger` from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded `BinaryLogger` object
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static BinaryLogger load(ARCHIVE& archive)
    {
        BinaryLogger logger;

        std::string output_directory;

        archive & output_directory
                & logger.cells_per_file
                & logger.cell_in_current_file
                & logger.current_file_number;

        logger.directory = output_directory;

        return logger;
    }

    /**
     * @brief The destructor
     */
    ~BinaryLogger();
};

}   // Simulation

}   // Drivers

}   // Races

#endif // __RACES_BINARY_LOGGER__
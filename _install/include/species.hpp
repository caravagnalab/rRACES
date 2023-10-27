/**
 * @file species.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines species representation
 * @version 0.18
 * @date 2023-10-13
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

#ifndef __RACES_SPECIES__
#define __RACES_SPECIES__

#include <iostream>

#include <map>
#include <vector>
#include <random>
#include <cstddef> // std::ptrdiff_t

#include "archive.hpp"
#include "time.hpp"
#include "cell.hpp"
#include "driver_genotype.hpp"


namespace Races 
{

namespace Drivers
{

namespace Simulation
{

class Tissue;

/**
 * @brief Cell species
 * 
 * This class represents the set of cells having the same driver genotype.
 */
class Species: public EpigeneticGenotype {
    /**
     * @brief A map from cell ids to the corresponding pointers
     */
    using CellIdToCell = std::map<CellId, CellInTissue *>;

    CellIdToCell cells;                 //!< species cells sorted by birth time/id
    CellIdToCell duplication_enabled;   //!< species cells that are duplication enabled

    Time last_insertion_time;     //!< the last insertion time

    /**
     * @brief A constant iterator for the species cells
     */
    class const_iterator {
        CellIdToCell::const_iterator it; //!< the cell pointer vector iterator

        const_iterator(const CellIdToCell::const_iterator it);
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   CellInTissue;
        using pointer           =   const CellInTissue*;
        using reference         =   const CellInTissue&;
        using iterator_category =   std::random_access_iterator_tag;

        const_iterator();

        inline reference operator*() const 
        { 
            return *(it->second); 
        }

        inline pointer operator->() 
        {
            return it->second;
        }

        // Prefix increment
        inline const_iterator& operator++() 
        {
            ++it;
            return *this;
        }  

        // Postfix increment
        const_iterator operator++(int);

        // Prefix decrement
        inline const_iterator& operator--() 
        {
            --it;
            return *this;
        }  

        // Postfix decrement
        const_iterator operator--(int);

        friend inline bool operator==(const const_iterator& a, const const_iterator& b)
        { 
            return a.it == b.it; 
        }

        friend inline bool operator!=(const const_iterator& a, const const_iterator& b)
        { 
            return a.it != b.it; 
        }

        friend class Species;
    };

    /**
     * @brief Compute the cell age from the last addition to the species
     * 
     * @param cell 
     * @return Time 
     */
    inline Time age(const CellInTissue *cell) const
    {
        return last_insertion_time-cell->get_birth_time();
    }

    /**
     * @brief Add a cell to the species
     * 
     * This method add the cell pointer by the parameter to the 
     * species. The species representation uses the memory pointed 
     * by the parameter.
     * 
     * @param cell is a pointer to the cell to be added to the species
     * @return pointer to the added cell
     */
    CellInTissue* add(CellInTissue* cell);

    /**
     * @brief Randomly select a cell among those in a id-cell pointer map
     * 
     * This method select a cell by using an exponential over 
     * cell ages.
     * 
     * @tparam GENERATOR is a random number generator type
     * @param generator is the random number generator used to select 
     *      the cell
     * @param id2cell_ptrs is the map id-cell pointers among which the 
     *      cell must be selected
     * @return a non-constant reference to a randomly selected
     *       cell in the species 
     */
    template<typename GENERATOR>
    const CellInTissue& choose_a_cell(GENERATOR& generator, 
                                      const CellIdToCell& id2cell_ptrs) const
    {
        if (id2cell_ptrs.size()==0) {
            throw std::domain_error("No cells can be selected");
        }

        std::exponential_distribution<Time> exp_dist;

        auto aimed_partial_ages_sum = exp_dist(generator);

        Time partial_ages_sum = 0;
        for (const auto& [cell_id, cell_ptr] : id2cell_ptrs) {
            if (partial_ages_sum == 0) {
                aimed_partial_ages_sum *= (10*age(cell_ptr));
            }

            partial_ages_sum += age(cell_ptr);

            if (partial_ages_sum>=aimed_partial_ages_sum) {
                return *cell_ptr;
            }
        }

        return *(cells.begin()->second);
    }

    /**
     * @brief An empty constructor
     */
    Species();

    /**
     * @brief Reset the species
     * 
     * This method resets the species and removes and deletes all 
     * its cells.
     */
    void reset();
public:
    /**
     * @brief A constructor
     * 
     * @param genotype is the epigenetic genotype of the species
     */
    explicit Species(const EpigeneticGenotype& genotype);

    /**
     * @brief A copy constructor
     * 
     * @param orig is the template object
     */
    Species(const Species& orig);

    /**
     * @brief The copy operator
     * 
     * @param orig is the original object 
     * @return a reference to the updated object
     */
    Species& operator=(const Species& orig);

    /**
     * @brief The replace operator
     * 
     * @param orig is the original object 
     * @return a reference to the updated object
     */
    Species& operator=(Species&& orig);

    /**
     * @brief Get the initial iterator for the species cells
     * 
     * @return the initial iterator for the species cells
     */
    const_iterator begin() const;

    /**
     * @brief Get the final iterator for the species cells
     * 
     * @return the final iterator for the species cells
     */
    const_iterator end() const;

    /**
     * @brief Get the number of cells in the species
     * 
     * @return the number of cells in the species
     */
    inline size_t num_of_cells() const
    {
        return cells.size();
    }

    /**
     * @brief Get the number of cells available for an event
     * 
     * @param event_type is the event for which the cells are needed
     * @return the number of cells available for `event_type`
     */
    size_t num_of_cells_available_for(const CellEventType& event_type) const;

    /**
     * @brief Randomly select a cell for a specific event type
     * 
     * This method select a cell by using an exponential over 
     * cell ages.
     * 
     * @tparam GENERATOR is a random number generator type
     * @param generator is the random number generator used to select 
     *      the cell
     * @param event_type is the type of the event that will occurs 
     *      to the selected cell
     * @return a non-constant reference to a randomly selected
     *       cell in the species 
     */
    template<typename GENERATOR>
    const CellInTissue& choose_a_cell(GENERATOR& generator,
                                      const CellEventType event_type) const
    {
        switch (event_type) {
            case CellEventType::DIE:
            case CellEventType::EPIGENETIC_EVENT:
            case CellEventType::DRIVER_GENETIC_MUTATION:
                return choose_a_cell(generator, cells);
            case CellEventType::DUPLICATE:
            case CellEventType::DUPLICATION_AND_EPIGENETIC_EVENT:
                return choose_a_cell(generator, duplication_enabled);
            default:
                throw std::domain_error("Unsupported event type");
        }
    }

    /**
     * @brief Randomly select a cell in the species
     * 
     * This method select a cell by using an exponential over 
     * cell ages.
     * 
     * @tparam GENERATOR is a random number generator type
     * @param generator is the random number generator used to select 
     *      the cell
     * @return a non-constant reference to a randomly selected
     *       cell in the species 
     */
    template<typename GENERATOR>
    inline const CellInTissue& choose_a_cell(GENERATOR& generator) const
    {
        return choose_a_cell(generator, cells);
    }

    /**
     * @brief Remove a cell from the species and delete it
     * 
     * @param cell_id is the id of the cell to be deleted
     */
    void erase(const CellId& cell_id);

    /**
     * @brief Add a cell to the species
     * 
     * By default the duplication of the cell is switched on.
     * 
     * @param cell is the cell to be added to the species
     * @return pointer to the added cell
     */
    CellInTissue* add(CellInTissue&& cell);

    /**
     * @brief Add a cell to the species
     * 
     * @param cell is the cell to be added to the species
     * @return pointer to the added cell
     */
    CellInTissue* add(CellInTissue& cell);

    /**
     * @brief Switch on/off the duplication of a cell
     * 
     * @param cell is the identifier of the cell on which 
     *          the duplication will be enabled/disabled
     * @param duplication_on is a Boolean flag to enable 
     *          (True)/disable(False) duplication on `cell` 
     */
    void switch_duplication_for(const CellId& cell_id,
                                const bool duplication_on);

    /**
     * @brief Enable the duplication of a cell
     * 
     * @param cell is the identifier of the cell to be enabled
     *          for the duplication
     */
    void enable_duplication_for(const CellId& cell_id);

    /**
     * @brief Disable the duplication of a cell
     * 
     * @param cell is the identifier of the cell to be disabled
     *          for the duplication
     */
    void disable_duplication_for(const CellId& cell_id);

    /**
     * @brief Get cell by identifier
     * 
     * @param cell_id is the identifier of the aimed cell
     * @return a non-constant reference to the aimed cell
     */
    CellInTissue& operator()(const CellId& cell_id);

    /**
     * @brief Get cell by identifier
     * 
     * @param cell_id is the identifier of the aimed cell
     * @return a constant reference to the aimed cell
     */
    const CellInTissue& operator()(const CellId& cell_id) const;

    /**
     * @brief Test whether a cell is already contained in a species
     * 
     * @param cell_id is the identifier of the cell to be tested
     * @return `true` is an only if the cell is already in the species
     */
    inline bool contains(const CellId& cell_id) const
    {
        return cells.find(cell_id) != cells.end();
    }

    /**
     * @brief The destroyer
     */
    ~Species();

    /**
     * @brief Save a species in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & static_cast<const EpigeneticGenotype &>(*this);

        archive & last_insertion_time;

        // save species cells
        archive & cells.size();
        for (const auto& [cell_id, cell_ptr]: cells) {
            archive & *cell_ptr;
        }

        // save duplicated ids
        archive & duplication_enabled.size();
        for (const auto& [cell_id, cell_ptr]: duplication_enabled) {
            archive & cell_id;
        }
    }

    /**
     * @brief Load a species from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded species
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static Species load(ARCHIVE& archive)
    {
        auto genotype = EpigeneticGenotype::load(archive);
        Species species(genotype);

        archive & species.last_insertion_time;

        size_t num_of_cells;

        // load species cells
        archive & num_of_cells;
        for (size_t i=0; i<num_of_cells; ++i) {
            CellInTissue *cell = new CellInTissue();

            archive & *cell;

            species.cells.insert({cell->get_id(), cell});
        }

        // load duplicated enabled ids
        archive & num_of_cells;
        for (size_t i=0; i<num_of_cells; ++i) {
            CellId cell_id;

            archive & cell_id;

            species.enable_duplication_for(cell_id);
        }

        return species;
    }

    friend class Simulation;

    friend void swap(Species& a, Species& b);
};


/**
 * @brief Swap two species
 * 
 * @param a is a species
 * @param b is a species
 */
void swap(Species& a, Species& b);

}   // Simulation

}   // Drivers

}   // Races

namespace std
{

/**
 * @brief Write the information about a species in an output stream
 * 
 * @param out is the output stream
 * @param species is the species to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Drivers::Species& species);

}   // std

#endif // __RACES_SPECIES__
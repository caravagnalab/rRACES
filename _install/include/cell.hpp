/**
 * @file cell.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines cell representation
 * @version 0.19
 * @date 2023-10-02
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

#ifndef __RACES_CELL__
#define __RACES_CELL__

#include <cstdint>

#include "archive.hpp"
#include "time.hpp"
#include "position.hpp"
#include "driver_genotype.hpp"

namespace Races
{

namespace Drivers
{

namespace Simulation
{

class Simulation;

}

/**
 * @brief The cell identifier type
 */
using CellId = uint64_t;

/**
 * @brief A class representing cells
 */
class Cell {
protected:
    static uint64_t counter;             //!< Total number of cell along the computation
    CellId id;                           //!< Cell identifier
    CellId parent;                       //!< Parent cell identifier

    Time birth_time;                     //!< Cell birth time

    EpigeneticGenotypeId epigenetic_id;  //!< cell species reference

    /**
     * @brief The empty constructor
     */
    Cell();

public:
    /**
     * @brief Create a new cell with no passenger mutations
     * 
     * @param epigenetic_id is the epigenetic genotype identifier
     */
    explicit Cell(const EpigeneticGenotypeId epigenetic_id);

    /**
     * @brief Create a new cell
     * 
     * @param epigenetic_id is the epigenetic genotype identifier
     * @param parent_id is the parent cell identifier
     */
    Cell(const EpigeneticGenotypeId epigenetic_id, const CellId parent_id);

    /**
     * @brief Create a new cell
     * 
     * @param epigenetic_id is the epigenetic genotype identifier
     * @param parent_id is the parent cell identifier
     * @param birth_time is the cell birth time
     */
    Cell(const EpigeneticGenotypeId epigenetic_id, const CellId parent_id, const Time birth_time);

    /**
     * @brief Get the cell identifier
     * 
     * @return a constant reference to the cell identifier
     */
    inline const CellId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the cell parent identifier
     * 
     * @return a constant reference to the cell parent identifier
     */
    inline const CellId& get_parent_id() const
    {
        return parent;
    }

    /**
     * @brief Get the cell birth time
     * 
     * @return the cell birth time 
     */
    inline const Time& get_birth_time() const
    {
        return birth_time;
    }

    /**
     * @brief Get the cell driver genotype
     * 
     * @return a constant reference to the cell driver genotype
     */
    inline const EpigeneticGenotypeId& get_epigenetic_id() const
    {
        return epigenetic_id;
    }

    /**
     * @brief Generate a descendent cell
     * 
     * @return the generated cell
     */
    inline Cell generate_descendent(const Time& time) const
    {
        return Cell(get_epigenetic_id(), get_id(), time);
    }

    /**
     * @brief Save a cell in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & id 
                & parent
                & birth_time
                & epigenetic_id;
    }

    /**
     * @brief Load a cell from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load cell
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static Cell load(ARCHIVE& archive)
    {
        Cell cell;

        archive & cell.id
                & cell.parent
                & cell.birth_time
                & cell.epigenetic_id;

        return cell;
    }

    friend class Tissue;
    friend class Species;
    friend class Simulation::Simulation; 

    friend void swap(Cell& a, Cell &b);
};

/**
 * @brief Swap two cells
 * 
 * @param a is a cell
 * @param b is a cell
 */
void swap(Cell& a, Cell &b);

/**
 * @brief Labelled cell
 * 
 * This template represents cells labelled by information, e.g., 
 * their birth time or a list of their passenger mutations. 
 * 
 * @tparam LABEL is the label type
 */
template<typename LABEL>
class LabelledCell : public Cell
{
    /**
     * @brief The empty constructor
     */
    LabelledCell():
        Cell(), label()
    {}
public:
    LABEL label;        //!< The cell label

    /**
     * @brief A constructor
     * 
     * @param cell is the unlabelled cell
     * @param label is the label
     */
    LabelledCell(const Cell& cell, LABEL label):
        Cell(cell), label(label)
    {}

    /**
     * @brief Save a labelled cell in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        static_cast<const Cell *>(this)->save(archive);

        archive & label;
    }

    /**
     * @brief Load a labelled cell from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the load labelled cell
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static LabelledCell<LABEL> load(ARCHIVE& archive)
    {
        LabelledCell<LABEL> cell;
        archive & static_cast<Cell &>(cell) 
                & cell.label;

        return cell;
    }

    template<typename LABEL2>
    friend void swap(LabelledCell<LABEL2>& a, LabelledCell<LABEL2> &b);
};

/**
 * @brief Swap two labelled cells
 * 
 * @param a is a cell
 * @param b is a cell
 */
template<typename LABEL>
void swap(LabelledCell<LABEL>& a, LabelledCell<LABEL> &b)
{
    std::swap(static_cast<Cell&>(a),static_cast<Cell&>(b));
    std::swap(a.label,b.label);
}

class Species;

namespace Simulation
{

/**
 * @brief A class to represent a cell in a tissue
 * 
 * A cell in a tissue is a cell provided with a position 
 * in the tissue
 */
class CellInTissue : public Cell, public PositionInTissue {

protected:
    /**
     * @brief The empty constructor
     */
    CellInTissue();

    /**
     * @brief A cell in tissue constructor
     * 
     * @param genotype is the cell driver genotype id
     * @param position is the cell position
     */
    CellInTissue(const EpigeneticGenotypeId genotype, const PositionInTissue& position);

    /**
     * @brief A cell in tissue constructor
     * 
     * @param cell is the cell
     * @param position is the cell position
     */
    CellInTissue(const Cell& cell, const PositionInTissue& position);

    /**
     * @brief A cell in tissue constructor
     * 
     * This constructor is meant to build non-driver-genotype cells.
     * 
     * @param position is the cell position
     */
    explicit CellInTissue(const PositionInTissue& position);

    /**
     * @brief A cell in tissue constructor
     * 
     * This constructor is meant to build non-driver-genotype cells.
     * 
     * @param x is the x axis position in the tissue
     * @param y is the y axis position in the tissue
     * @param z is the z axis position in the tissue
     */
    CellInTissue(const AxisPosition& x, const AxisPosition& y, const AxisPosition& z);

    /**
     * @brief Copy a cell into a cell in a tissue
     * 
     * @param position is the position in the tissue
     * @return a reference to the updated object
     */
    CellInTissue& operator=(const PositionInTissue& position);

public:

    /**
     * @brief Copy a cell into a cell in a tissue
     * 
     * @param cell is the original free cell
     * @return a reference to the updated object
     */
    CellInTissue& operator=(const Cell& cell);

    /**
     * @brief Test whether the cell has a driver genotype
     * 
     * @return true if and only if the cell has a driver genotype
     */
    inline bool has_driver_mutations() const
    {
        return get_epigenetic_id() != NON_DRIVER_GENOTYPE;
    }

    /**
     * @brief Save a cell in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & static_cast<const Cell &>(*this) & x & y & z;
    }

    /**
     * @brief Load a cell from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded cell
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static CellInTissue load(ARCHIVE& archive)
    {
        CellInTissue cell;

        archive & static_cast<Cell &>(cell) 
                & cell.x & cell.y & cell.z;

        return cell;
    }

    friend class Tissue;
    friend class Species;
    friend class CellInTissueProxy;
};

}   // Simulation

}   // Drivers

}   // Races

namespace std
{

/**
 * @brief Write a cell in an output stream
 * 
 * @param os is the output stream
 * @param cell is the cell to be streamed
 * @return a reference to the updated stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Drivers::Cell& cell);

/**
 * @brief Write a cell in a tissue in an output stream
 * 
 * @param os is the output stream
 * @param cell is the cell in the tissue to be streamed
 * @return a reference to the updated stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Drivers::Simulation::CellInTissue& cell);

} // std

#endif // __RACES_CELL__
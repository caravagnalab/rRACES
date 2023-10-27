/**
 * @file position_set.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to represent tissue position set
 * @version 0.5
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

#ifndef __RACES_POSITION_SET__
#define __RACES_POSITION_SET__

#include <cstdint>
#include <vector>

#include "tissue.hpp"
#include "archive.hpp"

namespace Races
{

namespace Drivers
{

/**
 * @brief The basic class for position set
 */
struct BasicPositionSet
{
    struct const_iterator {
    };

    /**
     * @brief Get the initial sampled cell iterator
     * 
     * This method must be over-loaded.
     * 
     * @return the initial sampled cell iterator
     * @throws std::runtime_error
     */
    const_iterator begin() const;

    /**
     * @brief Get the final sampled cell iterator
     * 
     * This method must be over-loaded.
     * 
     * @return the final sampled cell iterator
     * @throws std::runtime_error
     */
    const_iterator end() const;

    /**
     * @brief Return the number of positions in the rectangle
     * 
     * @return the number of positions in the rectangle
     */
    virtual size_t size() const = 0;

    /**
     * @brief The destroyer
     */
    virtual ~BasicPositionSet();
};

/**
 * @brief A hyper-rectangle of positions
 */
struct RectangleSet : public BasicPositionSet
{
    Simulation::PositionInTissue lower_corner;  //!< The lower corner in the hyper-rectangle
    Simulation::PositionInTissue upper_corner;  //!< The upper corner in the hyper-rectangle

    /**
     * @brief Constant iterators for the sample cells
     */
    class const_iterator
    {
        const RectangleSet* rectangle;

        Simulation::PositionInTissue pos;

        /**
         * @brief A constructor
         * 
         * @param rectangle is the a rectangle set
         * @param position is the current position of the new iterator
         */
        const_iterator(const RectangleSet* rectangle, const Simulation::PositionInTissue& position);
    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   Simulation::PositionInTissue;
        using pointer           =   const Simulation::PositionInTissue*;
        using reference         =   const Simulation::PositionInTissue&;
        using iterator_category =   std::bidirectional_iterator_tag;

        /**
         * @brief An empty constructor
         */
        const_iterator();

        /**
         * @brief Reference operator
         * 
         * @return a reference to the species pointer by the iterator 
         */
        inline reference operator*() const 
        { 
            return pos;
        }

        /**
         * @brief Pointer operator
         * 
         * @return a pointer to the species pointer by the iterator 
         */
        inline pointer operator->() 
        {
            return &(pos);
        }

        /**
         * @brief The prefix increment
         * 
         * @return a reference to the updated object
         */
        const_iterator& operator++();

        /**
         * @brief The postfix increment
         * 
         * @return a copy of the original object
         */
        inline const_iterator operator++(int)
        {
            RectangleSet::const_iterator copy(*this);

            this->operator++();

            return copy;
        }

        /**
         * @brief The prefix decrement
         * 
         * @return a reference to the updated object
         */
        const_iterator& operator--();

        /**
         * @brief The postfix decrement
         * 
         * @return a copy of the original object
         */
        inline const_iterator operator--(int)
        {
            RectangleSet::const_iterator copy(*this);

            this->operator--();

            return copy;
        }

        /**
         * @brief Test whether two iterators are the same
         * 
         * @param a is the first iterator to compare
         * @param b is the second iterator to compare
         * @return `true` if and only if the two iterators 
         *      refer to the same object
         */
        friend inline bool operator==(const const_iterator& a, const const_iterator& b)
        { 
            return (a.pos == b.pos) && (a.rectangle == b.rectangle); 
        }

        friend struct RectangleSet;
    };

    /**
     * @brief The empty constructor
     */
    RectangleSet();

    /**
     * @brief A punctual hyper-rectangle set constructor
     * 
     * @param position is the only position in the hyper-rectangle set
     */
    RectangleSet(const Simulation::PositionInTissue& position);

    /**
     * @brief A hyper-rectangle set constructor
     * 
     * @param lower_corner is the hyper-rectangle lower corner
     * @param upper_corner is the hyper-rectangle upper corner
     */
    RectangleSet(const Simulation::PositionInTissue& lower_corner, 
                 const Simulation::PositionInTissue& upper_corner);

    /**
     * @brief A cuboid set constructor
     * 
     * @param lower_corner is the hyper-rectangle lower corner 
     * @param x_size is the hyper-rectangle size along the x-axis
     * @param y_size is the hyper-rectangle size along the y-axis
     * @param z_size is the hyper-rectangle size along the z-axis
     */
    RectangleSet(const Simulation::PositionInTissue& lower_corner, 
                 const Simulation::AxisSize& x_size, const Simulation::AxisSize& y_size, 
                 const Simulation::AxisSize& z_size);

    /**
     * @brief A rectangle set constructor
     * 
     * @param lower_corner is the hyper-rectangle lower corner 
     * @param x_size is the hyper-rectangle size along the x-axis
     * @param y_size is the hyper-rectangle size along the y-axis
     */
    RectangleSet(const Simulation::PositionInTissue& lower_corner, 
                 const Simulation::AxisSize& x_size, const Simulation::AxisSize& y_size);

    /**
     * @brief Get the initial position iterator
     * 
     * @return the initial position iterator
     */
    const_iterator begin() const;

    /**
     * @brief Get the final position iterator
     * 
     * @return the final position iterator
     */
    const_iterator end() const;

    /**
     * @brief Return the number of positions in the rectangle
     * 
     * @return the number of positions in the rectangle
     */
    size_t size() const;

    /**
     * @brief Save a rectangle set
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & lower_corner 
                & upper_corner;
    }

    /**
     * @brief Load a timed genomic mutation from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded timed genomic mutation
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static RectangleSet load(ARCHIVE& archive)
    {
        Simulation::PositionInTissue lower_corner, upper_corner;

        archive & lower_corner
                & upper_corner;

        return {lower_corner, upper_corner};
    }

    /**
     * @brief Test whether two position rectangle sets are the same
     * 
     * @param a is the first position rectangle set to compare
     * @param b is the second position rectangle set to compare
     * @return `true` if and only if the two rectangle sets 
     *      refer to the same region
     */
    friend inline bool operator==(const RectangleSet& a, const RectangleSet& b)
    { 
        return ((a.lower_corner==b.lower_corner)
                && (a.upper_corner==b.upper_corner));
    }

    friend RectangleSet::const_iterator;
};

/**
 * @brief Test whether two iterators differs
 * 
 * @param a is the first iterator to compare
 * @param b is the second iterator to compare
 * @return `true` if and only if the two iterators 
 *      do not refer to the same object
 */
inline bool operator!=(const RectangleSet::const_iterator& a, const RectangleSet::const_iterator& b)
{ 
    return !(a==b); 
}

}   // Drivers

}   // Races

namespace std {
/**
 * @brief Write the data of a rectangle position set in a stream
 * 
 * @param os is the output stream
 * @param rectangle is a rectangle position set
 * @return a reference to the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Drivers::RectangleSet& rectangle);

}   // std

#endif // __RACES_POSITION_SET__
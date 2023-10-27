/**
 * @file phyloXML.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines phyloXML stream
 * @version 0.5
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

#ifndef __RACES_PHYLOXML__
#define __RACES_PHYLOXML__

#include <ostream>

#include "phylogenetic_forest.hpp"
#include "palette.hpp"

namespace Races
{

namespace Drivers
{

/**
 * @brief The namespace of IO classes in driver simulation
 */
namespace IO
{

/**
 * @brief A phyloXML stream for PhylogeneticForest
 */
class phyloXMLStream
{
    std::ostream& os;   //!< the output stream

    bool closed;        //!< a flag to signal whether the stream has already been closed

    size_t indent_level;            //!< the indentation level
    std::string indent_symbols;     //!< the indentation symbols

    std::string indent_string;      //!< the indentation string

    /**
     * @brief Stream the cell taxonomy
     * 
     * @param genotype_id is the cell genotype id
     * @return a reference to the updated object
     */
    phyloXMLStream& operator<<(const EpigeneticGenotypeId& genotype_id);

    /**
     * @brief Stream the branch color
     * 
     * @param color is the branch color
     * @return a reference to the updated object
     */
    phyloXMLStream& operator<<(const UI::Color& color);
public:
    /**
     * @brief Create a new phyloXML stream
     * 
     * @param os is the underline output stream
     * @param indentation_symbols are the indentation symbols
     */
    explicit phyloXMLStream(std::ostream& os=std::cout, const std::string& indentation_symbols=" ");

    /**
     * @brief  Stream a phylogenetic forest
     * 
     * @param forest is the forest to stream
     * @return a reference to the updated object
     */
    phyloXMLStream& operator<<(const PhylogeneticForest& forest);

    /**
     * @brief  Stream a phylogenetic forest node
     * 
     * @param node is the node to stream
     * @return a reference to the updated object
     */
    phyloXMLStream& operator<<(const PhylogeneticForest::const_node& node);

    /**
     * @brief Close the phyloXML tag
     * 
     * This method close the phyloxml tag and does not handle the 
     * underline basic stream.
     */
    void close();

    /**
     * @brief Change the indentation level
     * 
     * @param level is the new indentation level
     */
    void change_indentation_level(const size_t level);

    /**
     * @brief The destroyer
     */
    ~phyloXMLStream();
};

}   // IO

}   // Drivers

}   // Race

#endif // __RACES_PHYLOXML__
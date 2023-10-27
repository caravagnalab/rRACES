/**
 * @file fasta_utils.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines support utilities for FASTA files
 * @version 0.7
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

#ifndef __RACES_FASTA_UTILS__
#define __RACES_FASTA_UTILS__

#include <list>
#include <string>
#include <memory>

#include "genomic_position.hpp"

namespace Races
{

/**
 * @brief The generic Input/Output namespace
 * 
 */
namespace IO
{

namespace FASTA
{

/**
 * @brief A structure to decode chromosome sequences name
 */
struct SeqNameDecoder
{
    /**
     * @brief Establish whether a sequence is a chromosome by its name
     * 
     * @param seq_name is a FASTA sequence name
     * @param chr_id is the variable where the chromosome id will be placed
     * @return `true` if and only if `seq_name` correspond to a DNA chromosome sequence name
     */
    virtual bool is_chromosome_header(const std::string& seq_name, Passengers::ChromosomeId& chr_id) const = 0;

    /**
     * @brief The destroyer
     */
    virtual ~SeqNameDecoder();
};

/**
 * @brief A structure to decode Ensembl chromosome sequences name 
 */
struct EnsemblSeqNameDecoder : public SeqNameDecoder
{
    /**
     * @brief Establish whether a sequence is a chromosome by its name
     * 
     * @param seq_name is a FASTA sequence name
     * @param chr_id is the variable where the chromosome id will be placed
     * @return `true` if and only if `seq_name` correspond to a DNA chromosome sequence 
     *      name in Ensembl format
     */
    bool is_chromosome_header(const std::string& seq_name, Passengers::ChromosomeId& chr_id) const override;
};

/**
 * @brief A structure to decode NCBI chromosome sequences name 
 */
struct NCBISeqNameDecoder : public SeqNameDecoder
{
    /**
     * @brief Establish whether a sequence is a chromosome by its name
     * 
     * @param seq_name is a FASTA sequence name
     * @param chr_id is the variable where the chromosome id will be placed
     * @return `true` if and only if `seq_name` correspond to a DNA chromosome sequence 
     *      name in NCBI format
     */
    bool is_chromosome_header(const std::string& seq_name, Passengers::ChromosomeId& chr_id) const override;
};

/**
 * @brief The list of sequence name decorders
 * 
 * The elements in this list are used by the function
 * `Races::IO::FASTA::is_chromosome_header(...)`.
 */
extern std::list<std::shared_ptr<SeqNameDecoder>> seq_name_decoders;

/**
 * @brief Establish whether a sequence is a chromosome by its name
 * 
 * @param seq_name is a FASTA sequence name
 * @param chr_id is the variable where the chromosome id will be placed
 * @return `true` if and only if `seq_name` correspond to a DNA chromosome sequence 
 *      name in NCBI format
 */
bool is_chromosome_header(const std::string& seq_name, Passengers::ChromosomeId& chr_id);

}   // FASTA

}   // IO

}   // Races

#endif // __RACES_FASTA_UTILS__
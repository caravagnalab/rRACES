/**
 * @file cna.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class for copy number alterations
 * @version 0.3
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

#ifndef __RACES_CNA__
#define __RACES_CNA__

#include "allele.hpp"
#include "genomic_region.hpp"

namespace Races
{

namespace Passengers
{

/**
 * @brief Copy Number Alteration type
 */
struct CopyNumberAlteration
{
    /**
     * @brief The CNA type
     */
    enum class Type {
        AMPLIFICATION,
        DELETION
    };

    GenomicRegion region; //!< CNA region
    AlleleId source;      //!< source allele id
    AlleleId dest;        //!< destination allele id
    Type type;            //!< amplification/deletion flag

    /**
     * @brief The empty constructor
     */
    CopyNumberAlteration();

    /**
     * @brief A constructor
     * 
     * @param region is the region in which CNA occurred
     * @param source is the allele from which the region is copied/removed
     * @param destination is the allele in which the region is copied/removed
     * @param type is the CNA type
     */
    CopyNumberAlteration(const GenomicRegion& region, const AlleleId& source, const AlleleId& destination, const Type& type);

    /**
     * @brief Build a new amplification
     * 
     * @param region is the region in which CNA occurred
     * @param source is the allele from which the region is copied
     * @param destination is the allele in which the region is copied
     */
    static inline CopyNumberAlteration new_amplification(const GenomicRegion& region,
                                                          const AlleleId& source,
                                                          const AlleleId& destination)
    {
        return CopyNumberAlteration(region, source, destination, Type::AMPLIFICATION);
    }

    /**
     * @brief Build a new amplification
     * 
     * @param region is the region in which CNA occurred
     * @param allele is the allele from which the region is removed
     */
    static inline CopyNumberAlteration new_deletion(const GenomicRegion& region, const AlleleId& allele)
    {
        return CopyNumberAlteration(region, allele, allele, Type::DELETION);
    }
};

}   // Passengers

}   // Races

#endif // __RACES_CNA__
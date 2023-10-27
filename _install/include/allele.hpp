/**
 * @file allele.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines allele representation
 * @version 0.6
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

#ifndef __RACES_ALLELE__
#define __RACES_ALLELE__

#include <map>

#include "snv.hpp"
#include "genomic_region.hpp"

namespace Races
{

namespace Passengers
{

/**
 * @brief A identifier type for alleles
 */
typedef size_t AlleleId;

/**
 * @brief A class to represent a fragment of an allele
 * 
 */
class AlleleFragment : public GenomicRegion
{
    std::map<GenomicPosition, SNV> snvs;  //!< the fragment SNVs

    /**
     * @brief Split an allele fragment
     * 
     * This method cuts an allele fragment in a position, 
     * generates a new allele fragment beginning at the specified 
     * position, and updates the length of the original one so 
     * that the two fragments are contiguous, i.e., one of the two
     * follows the other one. 
     * The SNVs of the original allele fragment are distributed on
     * the two allele fragments according to their genomic position.
     * 
     * @param split_point is the position of the new allele fragment
     * @return the allele fragment originated by the split
     * @throw std::domain_error the allele fragment does not contain 
     *          `split_point` or `split_point` and the allele 
     *          fragment initial point are the same
     */
    AlleleFragment split(const GenomicPosition& split_point);

    /**
     * @brief Split an allele fragment
     * 
     * This method cuts an allele fragment in a position, 
     * generates a new allele fragment beginning at the specified 
     * position, and updates the length of the original one so 
     * that the two fragments are contiguous, i.e., one of the two
     * follows the other one. 
     * The SNVs of the original allele fragment are distributed on
     * the two allele fragments according to their genomic position.
     * 
     * @param split_point is the position of the new allele fragment
     * @return the allele fragment originated by the split
     * @throw std::domain_error the allele fragment does not contain 
     *          `split_point` or `split_point` and the allele 
     *          fragment initial point are the same
     */
    inline AlleleFragment split(GenomicPosition&& split_point)
    {
        return split(split_point);
    }

public:
    /**
     * @brief The allele fragment length
     */
    using Length = GenomicRegion::Length;

    /**
     * @brief The empty constructor
     */
    AlleleFragment();

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the chromosome identifier
     * @param begin is the initial position of the allele fragment
     * @param end is the final position of the allele fragment
     */
    AlleleFragment(const ChromosomeId& chromosome_id,
                   const ChrPosition& begin, const ChrPosition& end);

    /**
     * @brief A constructor
     * 
     * @param genomic_region is the genomic region of the allele fragment
     */
    explicit AlleleFragment(const GenomicRegion& genomic_region);

    /**
     * @brief Get the allele fragment SNVs
     * 
     * @return a constant reference to the allele fragment SNVs
     */
    inline const std::map<GenomicPosition, SNV>& get_SNVs() const
    {
        return snvs;
    }

    /**
     * @brief Check whether the fragment context is free
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if the allele fragment does not 
     *      contains SNVs in the one-base neighborhood of 
     *      `genomic_position`
     */
    bool has_context_free(const GenomicPosition& genomic_position) const;

    /**
     * @brief Insert a new SNV
     * 
     * This method tries to insert a new SNV. It succeeds 
     * if no other SNVs are contained in the context.
     * 
     * @param snv is the SNV to insert
     * @return `true` if and only if the SNV insertion 
     *          has succeeded
     */
    bool insert(const SNV& snv);

    /**
     * @brief Remove a SNV 
     * 
     * This method tries to remove a SNV from a genomic position. 
     * It succeeds if the allele fragment contains a SNV in the 
     * specified position.
     * 
     * @param genomic_position is the genomic position containing the 
     *          SVN to be removed
     * @return `true` if and only if the removal has succeeded
     */
    bool remove_SNV(const GenomicPosition& genomic_position);

    /**
     * @brief Copy part of an allele fragment
     * 
     * @param genomic_region is the genomic region to copy
     * @return an allelic fragment corresponding to `genomic_region` 
     *      that contains all the original allele fragment SNVs 
     *      laying in `genomic_region`
     */
    AlleleFragment copy(const GenomicRegion& genomic_region) const;

    friend class Allele;
};

/**
 * @brief A possible representation for an allele
 */
class Allele
{
    std::map<GenomicPosition, AlleleFragment> fragments;    //!< the sequence fragments
public:
    /**
     * @brief The allele length
     */
    using Length = AlleleFragment::Length;

    /**
     * @brief The empty constructor
     */
    Allele();

    /**
     * @brief A constructor
     * 
     * @param chromosome_id is the chromosome identifier
     * @param begin is the initial position of the allele
     * @param end is the final position of the allele
     */
    Allele(const ChromosomeId& chromosome_id, const ChrPosition& begin, const ChrPosition& end);

    /**
     * @brief A constructor
     * 
     * @param genomic_region is the genomic region of the allele
     */
    explicit Allele(const GenomicRegion& genomic_region);

    /**
     * @brief Get the allele fragments
     * 
     * @return a constant reference to the allele fragments 
     */
    inline const std::map<GenomicPosition, AlleleFragment>& get_fragments() const
    {
        return fragments;
    }

    /**
     * @brief Test whether the allele contains a genomic position
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if one of the allele fragments 
     *          contains `genomic_position`
     */
    bool contains(const GenomicPosition& genomic_position) const;

    /**
     * @brief Test whether the allele contains a genomic region
     * 
     * @param genomic_region is a genomic region
     * @return `true` if and only if one of the allele fragments 
     *          contains `genomic_regions`
     */
    bool contains(const GenomicRegion& genomic_region) const;

    /**
     * @brief Check whether the fragment context is free
     * 
     * @param genomic_position is a genomic position
     * @return `true` if and only if the allele fragment does not 
     *      contains SNVs in the one-base neighborhood of 
     *      `genomic_position`
     */
    bool has_context_free(const GenomicPosition& genomic_position) const;

    /**
     * @brief Insert a new SNV
     * 
     * This method tries to insert a new SNV. It succeeds 
     * if no other SNVs are contained in the context.
     * 
     * @param snv is the SNV to insert
     * @return `true` if and only if the SNV insertion 
     *          has succeeded
     */
    bool insert(const SNV& snv);

    /**
     * @brief Remove a SNV 
     * 
     * This method tries to remove a SNV from a genomic position. 
     * It succeeds if the allele fragment contains a SNV in the 
     * specified position.
     * 
     * @param genomic_position is the genomic position containing the 
     *          SVN to be removed
     * @return `true` if and only if the removal has succeeded
     */
    bool remove_SNV(const GenomicPosition& genomic_position);

    /**
     * @brief Copy part of an allele
     * 
     * @param genomic_region is the genomic region to copy
     * @return an allelic fragment corresponding to `genomic_region` 
     *      that contains all the original allele SNVs 
     *      laying in `genomic_region`
     */
    Allele copy(const GenomicRegion& genomic_region) const;

    /**
     * @brief Remove part of an allele
     * 
     * This method tries to remove the part of an allele 
     * corresponding to a genomic region. It succeeds only if the 
     * original allele contains the specified genomic 
     * region.
     * 
     * @param genomic_region is the genomic region to copy
     * @return `true` if and only if the allele contains 
     *          `genomic_region`
     */
    bool remove(const GenomicRegion& genomic_region);

    /**
     * @brief Get the allele SNVs
     * 
     * @return the allele SNVs
     */
    std::map<GenomicPosition, SNV> get_SNVs() const;

    /**
     * @brief Get the size of the allele
     * 
     * @return the size of the allele 
     */
    Length size() const;
};

}   // Passengers

}   // Races

namespace std
{

/**
 * @brief Write allele fragment data in a stream
 * 
 * @param os is the output stream
 * @param allele_fragment is the allele fragment to be written
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Passengers::AlleleFragment& allele_fragment);


/**
 * @brief Write allele data in a stream
 * 
 * @param os is the output stream
 * @param allele is the allele to be written
 * @return a reference to output stream
 */
std::ostream& operator<<(std::ostream& os, const Races::Passengers::Allele& allele);

} // std

#endif // __RACES_ALLELE__
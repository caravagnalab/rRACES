/**
 * @file phylogenetic_forest.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for phylogenetic forests
 * @version 0.13
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

#ifndef __RACES_PHYLOGENETIC_FOREST__
#define __RACES_PHYLOGENETIC_FOREST__

#include <map>
#include <set>
#include <vector>
#include <queue>

#include "cell.hpp"
#include "binary_logger.hpp"

namespace Races
{

namespace Drivers 
{

/**
 * @brief A class representing phylogenetic forests
 */
class PhylogeneticForest
{
    std::map<CellId,  Cell> cells;                  //!< The forest cell id-cell map
 
    std::set<CellId> roots;                         //!< The cell ids of the forest roots

    std::map<CellId, std::set<CellId>> branches;    //!< The descendant branches

    std::map<EpigeneticGenotypeId, GenotypeId> genotype_map;    //!< A map associating every epigenetic genotype to the corresponding genotype id 

    std::vector<Simulation::TissueSample> samples;  //!< The vector of the samples that produced the forest
    std::map<CellId, Simulation::TissueSample const *> coming_from;       //!< The map associating each leaf to sample which it comes from

    /**
     * @brief Grow a forest from a sample of cells
     * 
     * This function grows a forest from a sample of cells.
     * The sample cells are the leaves of the forest and 
     * their ancestors are loaded from a cell storage.
     * 
     * @tparam CELL_STORAGE is the type of the cell storage
     * @param sample_ids is the cell id sample
     * @param cell_storage is the cell storage
     * @param genotypes is the vector of epigenetic genotypes
     */
    template<typename CELL_STORAGE>
    void grow_forest_from(const std::list<CellId>& sample, 
                          CELL_STORAGE& cell_storage,
                          const std::vector<EpigeneticGenotype>& genotypes)
    {
        clear();

        for (const auto& epigenetic_genotype: genotypes) {
            genotype_map[epigenetic_genotype.get_id()] = epigenetic_genotype.get_genomic_id();
        }

        std::set<CellId> parent_ids;
        for (const auto& cell_id: sample) {
            parent_ids.insert(cell_id);

            // record leaves children, i.e., none
            branches[cell_id] = std::set<CellId>();
        }

        std::priority_queue<CellId> queue(parent_ids.begin(), parent_ids.end());

        while (!queue.empty()) {
            auto cell = cell_storage[queue.top()];

            queue.pop();
            parent_ids.erase(cell.get_id());

            // if the cell is not an initial cell
            if (cell.get_id()!=cell.get_parent_id()) {

                // the cell id is not in the queue
                if (parent_ids.count(cell.get_parent_id())==0) {

                    // add its id to the queue
                    queue.push(cell.get_parent_id());
                    parent_ids.insert(cell.get_parent_id());
                }
                branches[cell.get_parent_id()].insert(cell.get_id());
            } else {
                // it is a root
                roots.insert(cell.get_id());
            }

            cells.insert(std::make_pair(cell.get_id(),cell));
        }
    }
public:

    /**
     * @brief A constant node of the forest
     */
    class const_node 
    {
        PhylogeneticForest* forest;         //!< A pointer to the forest
        CellId cell_id;                     //!< The cell id of a cell in the forest

        /**
         * @brief A constructor for a constant node
         * 
         * @param forest is the forest of the node
         * @param cell_id is the cell id of a cell in the forest
         */
        const_node(const PhylogeneticForest* forest, const CellId cell_id);
    public:
        /**
         * @brief Cast to Cell
         * 
         * This method returns a constant reference to a cell. 
         * Notice that the cells in the forest should be modified 
         * exclusively by using `PhylogeneticForest::node` methods
         * 
         * @return a constant reference to a cell
         */
        inline operator const Cell&() const
        {
            return forest->cells.at(cell_id);
        }

        /**
         * @brief Get the parent node
         * 
         * @return the parent node of this node
         */
        const_node parent() const;

        /**
         * @brief Get the node direct descendants
         * 
         * @return a vector containing the direct descendants of
         *         this node
         */
        std::vector<const_node> children() const;

        /**
         * @brief Test whether this node is a leaf
         * 
         * @return `true` if and only if this node is a leaf
         */
        inline bool is_leaf() const
        {
            return forest->branches.at(cell_id).size()==0;
        }

        /**
         * @brief Test whether this node is a root
         * 
         * @return `true` if and only if this node is a root
         */
        inline bool is_root() const
        {
            const auto& cell = forest->cells.at(cell_id);

            return cell.get_id() == cell.get_parent_id();
        }

        /**
         * @brief Get the cell id of the node
         * 
         * @return the cell id of the node 
         */
        inline CellId get_id() const
        {
            return cell_id;
        }

        /**
         * @brief Get the sample that collected the cell 
         * 
         * @return the sample that collected the identifier 
         *      of this node in the tissue
         * @throw `std::domain_error` when the node is not 
         *      a leaf
         */
        const Simulation::TissueSample& get_sample() const;

        /**
         * @brief Get the node epigenetic genotype id
         * 
         * @return the node epigenetic genotype id
         */
        inline EpigeneticGenotypeId get_epigenetic_id() const
        {
            return forest->cells.at(cell_id).get_epigenetic_id();
        }

        /**
         * @brief Get the node genotype id
         * 
         * @return the node genotype id
         */
        GenotypeId get_genotype_id() const;

        /**
         * @brief Get the node forest
         * 
         * @return a constant reference to the node forest
         */
        inline const PhylogeneticForest& get_forest() const
        {
            return *forest;
        }

        friend class PhylogeneticForest;
    };

    /**
     * @brief A non-constant node of the forest
     */
    class node : public const_node
    {

        node(PhylogeneticForest* forest, const CellId cell_id);
    public:
        /**
         * @brief Cast to Cell
         * 
         * @return a non-constant reference to a cell
         */
        inline operator Cell&()
        {
            return forest->cells.at(cell_id);
        }

        /**
         * @brief Get the parent node
         * 
         * @return the parent node of this node
         */
        node parent();

        /**
         * @brief Get the node direct descendants
         * 
         * @return a vector containing the direct descendants of
         *         this node
         */
        std::vector<node> children();

        friend class PhylogeneticForest;
    };

    /**
     * @brief The empty constructor
     */
    PhylogeneticForest();

    /**
     * @brief Construct a phylogenetic forest by using a tissue sample of a simulation
     * 
     * This method builds a phylogenetic forest by using clone simulation 
     * pre-sampled cells as leaves.
     * 
     * @param simulation is a simulation
     */
    PhylogeneticForest(const Simulation::Simulation& simulation);

    /**
     * @brief Construct a phylogenetic forest by using a tissue sample of a simulation
     * 
     * @param simulation is a simulation
     * @param tissue_samples is a list of tissue samples coming from the simulation
     */
    PhylogeneticForest(const Simulation::Simulation& simulation,
                       const std::list<Simulation::TissueSample>& tissue_samples);

    /**
     * @brief Get a constant node with the given id
     * 
     * @param cell_id is the id of the aimed cell node
     * @return the corresponding constant cell node
     */
    const_node get_node(const CellId& cell_id) const;

    /**
     * @brief Get a node with the given id
     * 
     * @param cell_id is the id of the aimed cell node
     * @return the corresponding cell node
     */
    node get_node(const CellId& cell_id);

    /**
     * @brief Get the forest roots
     * 
     * @return std::vector<const_node> 
     */
    std::vector<const_node> get_roots() const;

    /**
     * @brief Get the forest roots
     * 
     * @return std::vector<const_node> 
     */
    std::vector<node> get_roots();

    /**
     * @brief Get the number of forest nodes
     * 
     * @return the number of forest nodes 
     */
    inline size_t num_of_nodes() const
    {
        return cells.size();
    }

    /**
     * @brief Get the tissue samples that produced the forest
     * 
     * @return a constant reference to a vector of tissue samples 
     *      that produced the phylogenetic forest
     */
    inline const std::vector<Simulation::TissueSample>& get_samples() const
    {
        return samples;
    }

    /**
     * @brief Clear the forest
     */
    void clear();
};

}   // Drivers

}   // Races

#endif // __RACES_PHYLOGENETIC_TREE__
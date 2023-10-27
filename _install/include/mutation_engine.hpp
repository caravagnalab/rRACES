/**
 * @file mutation_engine.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to place passenger mutations on the nodes of a phylogenetic forest
 * @version 0.14
 * @date 2023-10-24
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

#ifndef __RACES_MUTATION_ENGINE__
#define __RACES_MUTATION_ENGINE__

#include <map>
#include <stack>
#include <random>
#include <ostream>

#include "phylogenetic_forest.hpp"
#include "driver_genotype_id.hpp"

#include "context_index.hpp"
#include "genome_mutations.hpp"
#include "snv_signature.hpp"
#include "mutational_properties.hpp"

#include "filter.hpp"


namespace Races
{

namespace Passengers
{

/**
 * @brief Mutation Statistics
 */
class MutationStatistics
{
    /**
     * @brief SNV statistics
     */
    struct SNVStatistics
    {
        size_t mutated_alleles;            //!< Number of mutated alleles
        size_t num_of_cells;               //!< Total number of cells containing the SNV

        /**
         * @brief The empty constructor
         */
        SNVStatistics();
    };

    struct SampleMutationStatistics
    {
        size_t num_of_cells;          //!< Number of recorded cells

        std::map<SNV, SNVStatistics> SNVs;      //!< SNVs
        std::list<CopyNumberAlteration> CNAs;   //!< CNAs

        /**
         * @brief The empty constructor
         */
        SampleMutationStatistics();
    };

    SampleMutationStatistics overall_statistics;    //!< Overall statistics
    std::map<std::string, SampleMutationStatistics> sample_statistics;  //!< Statistics per sample

public:
    /**
     * @brief The empty constructor
     */
    MutationStatistics();

    /**
     * @brief Record the SVNs of a cell genome
     * 
     * @param sample_name is the name of the sample from which the cell has 
     *          been extracted
     * @param cell_mutations are the mutations of the cell whose statistics
     *          is about to be recorded
     * @return a reference to the updated object
     */
    MutationStatistics& record(const std::string& sample_name,
                               const CellGenomeMutations& cell_mutations);

    /**
     * @brief Record the SVNs in a list of sample genomic mutations
     * 
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] progress_bar is a progress bar pointer
     * @return a reference to the updated object
     */
    MutationStatistics& 
    record(const std::list<Races::Passengers::SampleGenomeMutations>& mutations_list,
           UI::ProgressBar* progress_bar=nullptr);

    /**
     * @brief Record the SVNs in a list of sample genomic mutations
     * 
     * @param[in] mutations_list is a list of sample mutations
     * @param[in,out] progress_bar is a progress bar pointer
     * @return a reference to the updated object
     */
    inline MutationStatistics& 
    record(const std::list<Races::Passengers::SampleGenomeMutations>& mutations_list,
           UI::ProgressBar& progress_bar)
    {
        return record(mutations_list, &progress_bar);
    }

    /**
     * @brief Write in a stream a table summarizing SNV statistics
     * 
     * @param os is the output stream
     * @param separator is the column separator
     * @return a reference to the updated stream
     */
    std::ostream& write_SNVs_table(std::ostream& os, const char separator='\t');

    /**
     * @brief Write in a stream a table summarizing CNA statistics
     * 
     * @param os is the output stream
     * @param separator is the column separator
     * @return a reference to the updated stream
     */
    std::ostream& write_CNAs_table(std::ostream& os, const char separator='\t');
};

/**
 * @brief The mutation engine
 * 
 * The objects of this class place mutations on the cell genome according 
 * to a phylogenetic tree
 * 
 * @tparam GENOME_WIDE_POSITION is the type used to represent genome-wise position
 * @tparam RANDOM_GENERATOR is the type of random generator
 */
template<typename GENOME_WIDE_POSITION = uint32_t, typename RANDOM_GENERATOR = std::mt19937_64>
class MutationEngine
{
    /**
     * @brief The SBS inverse cumulative distributions
     * 
     * SBS's are probability distributions for mutational types. If we arbitrary assign
     * an order to mutational type, given a SBS and a mutational type T, we can compute
     * the probability of randomly choosing a mutational type smaller or equal to T. 
     * The inverse cumulative SBS type associated such a probability to T itself for 
     * all the mutational type T.
     */
    using InverseCumulativeSBS = std::map<double, Races::Passengers::MutationalType>;

    /**
     * @brief Coefficients for mutational signatures
     */
    using MutationalCoefficients = std::map<std::string, double>;

    /**
     * @brief The stack of contexts removed from a context index
     */
    using ContextStack = std::stack<std::pair<MutationalContext, GENOME_WIDE_POSITION>>;

    RANDOM_GENERATOR generator; //!< the random generator

    ContextIndex<GENOME_WIDE_POSITION> context_index;  //!< the genome context index
    ContextStack context_stack;                             //!< the stack of the contexts removed from `context_index`

    size_t num_of_alleles;                                  //!< number of alleles 

    std::map<Time, MutationalCoefficients> timed_mutational_coefficients;   //!< the timed mutational coefficients
    std::map<std::string, InverseCumulativeSBS> inv_cumulative_SBSs;  //!< the inverse cumulative SBS

    SpeciesMutationalProperties mutational_properties;  //!< the species mutational properties

    /**
     * @brief Select a random value in a set
     * 
     * @tparam T is the type of objects in the set
     * @param values is the set of values from which a random value must be selected
     * @return a random value among those in `values`
     */
    template<typename T>
    T select_random_value(const std::set<T>& values)
    {   
        std::uniform_int_distribution<> u_dist(0,values.size()-1);
        size_t pos = u_dist(generator);

        size_t index=0;
        auto it = values.begin();
        while (index!=pos) {
            ++it;
            ++index;
        }

        return *it;
    }

    /**
     * @brief Compute the inverse cumulative SBS
     
     * SBS's are probability distributions for mutational types. If we arbitrary assign
     * an order to mutational type, given a SBS and a mutational type T, we can compute
     * the probability of randomly choosing a mutational type smaller or equal to T. 
     * The inverse cumulative SBS type associated such a probability to T itself for 
     * all the mutational type T.
     * 
     * @param m_signature is a mutational signature (SBS)
     * @return the corresponding mutational signature
     */
    static InverseCumulativeSBS get_inv_cumulative_SBS(const MutationalSignature& m_signature)
    {
        InverseCumulativeSBS inv_cumulative_SBS;
        double cumulative_prob = 0;
        for (const auto& [m_type, prob]: m_signature) {
            if (prob != 0) {
                cumulative_prob += prob;

                inv_cumulative_SBS.insert({cumulative_prob, m_type});
            }
        }

        return inv_cumulative_SBS;
    }

    /**
     * @brief Get the mutational coefficients according to a cell birth time
     * 
     * The active mutational coefficients depend on the simulation time. They 
     * must be selected in agreement with cell birth times.
     * 
     * @param cell is the cell whose birth time select the mutational coefficients
     * @return a constant reference to mutational coefficients that are associated 
     *         to `cell` birth time
     */
    const MutationalCoefficients&
    get_active_mutational_coefficients(const Drivers::Cell& cell) const
    {
        auto mc_it = timed_mutational_coefficients.upper_bound(cell.get_birth_time());

        return (--mc_it)->second;
    }

    /**
     * @brief Select a random SNV 
     * 
     * This method selects a random SNV among whose having a specified 
     * mutational type and available in a context index. The selected 
     * SNV is extracted from the context index and inserted into a stack 
     * to revert the selection.
     * 
     * @param[in] m_type is the mutational type of the SNV to be selected
     * @param[in] cause is the SNV cause
     * @return a SNV whose type is `m_type` and which was available in 
     *          `context_index`
     */
    SNV select_SVN(const MutationalType& m_type, const std::string& cause)
    {
        using namespace Races::Passengers;

        MutationalContext context = m_type.get_context();
        MutationalContext complement_context = context.get_complement();

        size_t total_pos = context_index[context].size();
            total_pos += context_index[complement_context].size();

        size_t index;
        {
            std::uniform_int_distribution<> u_dist(0,total_pos-1);
            index = u_dist(generator);
        }

        char mutated_base = m_type.get_replace_base();
        if (index >= context_index[context].size()) {
            index -= context_index[context].size();
            context = complement_context;

            mutated_base = MutationalContext::get_complement(mutated_base);
        }

        GENOME_WIDE_POSITION pos = context_index.extract(context, index);

        context_stack.push({context, pos});

        auto genomic_pos = context_index.get_genomic_position(pos);

        return {genomic_pos, context, mutated_base, cause};
    }

    /**
     * @brief Randomly select the number of SNVs according to a Poisson distribution
     * 
     * @param genome_size is the genome size
     * @param mean_passenger_mutations is the mean passenger mutations per duplication
     * @return is the randomly selected number of SNVs
     */
    size_t number_of_SNVs(GenomeMutations::Length genome_size,
                          const double& mean_passenger_mutations)
    {
        auto SNVs_mean = static_cast<int>(genome_size*mean_passenger_mutations);

        std::poisson_distribution<> p_dist(SNVs_mean);

        return p_dist(generator);
    }

    /**
     * @brief Try to place a SNV
     * 
     * @param snv is the SNV to place
     * @param cell_mutations are the cell mutations
     * @return `true` if and only if the SNV placement has succeed. This 
     *          method returns `false` when the SNV cannot be placed 
     *          because its context is not free or no allele are available
     *          for it.
     */
    bool place_SNV(const SNV& snv, GenomeMutations& cell_mutations)
    {
        auto allele_having_pos = cell_mutations.get_alleles_containing(snv);
        if (allele_having_pos.size()>0) {
            AlleleId allele_id = select_random_value(allele_having_pos);

            // try to apply the selected SNV
            if (cell_mutations.insert(snv, allele_id)) {
                return true;
            }
        }

        return false;
    }

    /**
     * @brief Place the mutations associated to a cell in the phylogenetic tree
     * 
     * @tparam GENOME_WIDE_POSITION is the type used to represent genome-wise position
     * @param node is a phylogenetic tree node representing a cell
     * @param cell_mutations are the cell mutations
     */
    void place_SNVs(const Drivers::PhylogeneticForest::const_node& node,
                    GenomeMutations& cell_mutations)
    {
        auto num_of_SNVs = number_of_SNVs(cell_mutations.allelic_size(),
                                          mutational_properties.at(node.get_epigenetic_id()).mu);

        // get the active SBS coefficients
        const auto& mc = get_active_mutational_coefficients(node);

        // for each active SBS coefficient 
        for (const auto& [SBS_name, coefficient]: mc) {

            // get the inverse cumulative SBS
            const auto& inv_cumulative_SBS = inv_cumulative_SBSs.at(SBS_name);

            // evaluate how many of the SNVs occurred to the considered cells
            // are due to the current SBS
            size_t SBS_num_of_SNVs = coefficient*num_of_SNVs;

            std::uniform_real_distribution<> u_dist(0,1);

            // while some of the requested SNVs have not been set
            unsigned int new_SNVs = 0;
            while (new_SNVs < SBS_num_of_SNVs) {
                // randomly pick a mutational type according to the current SBS
                auto it = inv_cumulative_SBS.lower_bound(u_dist(generator));

                if (it != inv_cumulative_SBS.end()) {
                    // select a SNV among those having the picked mutational type 
                    // and that available in the context index 
                    SNV snv = select_SVN(it->second, SBS_name);

                    if (place_SNV(snv, cell_mutations)) {
                        // if the SNV has been successfully applied, 
                        // then increase the number of new SNVs
                        ++new_SNVs;
                    }
                }
            }
        }
    }

    /**
     * @brief Try to place a CNA
     * 
     * @param CNA is the CNA to place
     * @param cell_mutations are the cell mutations
     * @return `true` if and only if the CNA placement has succeed. This 
     *          method returns `false` when the CNA cannot be placed 
     *          because no allele are available for it.
     */
    bool place_CNA(const CopyNumberAlteration& CNA, GenomeMutations& cell_mutations)
    {
        switch(CNA.type) {
            case CopyNumberAlteration::Type::AMPLIFICATION:
                return cell_mutations.amplify_region(CNA.region, CNA.source);
            case CopyNumberAlteration::Type::DELETION:
                return cell_mutations.remove_region(CNA.region, CNA.source);
            default:
                throw std::runtime_error("Unsupported CNA type");
        }
    }

    /**
     * @brief Place the driver specific SNVs
     * 
     * @param node is a phylogenetic tree node representing a cell
     * @param mutations are the cell mutations
     * @param cell_statistics are the statistics of the cell mutations
     */
    void place_driver_specific_mutations(const Drivers::PhylogeneticForest::const_node& node,
                                         GenomeMutations& cell_mutations)
    {
        if (node.is_root() || node.get_genotype_id()!=node.parent().get_genotype_id()) {
            const auto& driver_mp = mutational_properties.at(node.get_epigenetic_id());

            for (auto snv : driver_mp.SNVs) {
                snv.cause = driver_mp.name;

                place_SNV(snv, cell_mutations);
            }

            for (auto CNA : driver_mp.CNAs) {
                place_CNA(CNA, cell_mutations);
            }
        }
    }

    /**
     * @brief Place the mutations on the genomes of a phylogenetic forest node
     * 
     * This method recursively places the mutations on the nodes of a phylogenetic
     * forest. All the genome mutations associated to the forest leaves are 
     * collected and saved in a container that partition them according 
     * to the cell sample. 
     * 
     * @param[in,out] sample_mutation_map is a map associating tissue sample ids to 
     *              the corresponding sample genomic mutations
     * @param[in] node is a phylogenetic tree node representing a cell
     * @param[in] ancestor_mutations is the genomic mutation of the ancestor
     * @param[in,out] visited_nodes is the number of visited nodes
     * @param[in,out] progress_bar is a progress bar pointer
     */
    void place_mutations(std::map<Drivers::Simulation::TissueSampleId, 
                                  SampleGenomeMutations*>& sample_mutation_map,
                         const Drivers::PhylogeneticForest::const_node& node,
                         const GenomeMutations& ancestor_mutations,
                         size_t& visited_nodes, UI::ProgressBar *progress_bar)
    {
        using namespace Races::Passengers;

        size_t context_stack_size = context_stack.size();
        
        CellGenomeMutations cell_mutations(static_cast<const Drivers::Cell&>(node), ancestor_mutations);

        place_driver_specific_mutations(node, cell_mutations);
        place_SNVs(node, cell_mutations);
        //place_CNAs(node, mutations);

        ++visited_nodes;
        if (progress_bar != nullptr) {

            size_t percentage = (100*visited_nodes)/(node.get_forest()).num_of_nodes();
            if (percentage>progress_bar->get_progress()) {
                progress_bar->set_progress(percentage);
            }
        }

        if (node.is_leaf()) {
            auto& node_sample = node.get_sample();
            auto& sample_genomic_mutations = *(sample_mutation_map.at(node_sample.get_id()));
            sample_genomic_mutations.mutations.push_back(cell_mutations);
        } else {
            for (const auto child: node.children()) {
                place_mutations(sample_mutation_map, child, cell_mutations,
                                visited_nodes, progress_bar);
            }
        }

        // reverse context index extractions
        while (context_stack.size() > context_stack_size) {
            reinsert_last_context();
        }
    }

    /**
     * @brief Remove the context stack top and insert it in the context index
     */
    void reinsert_last_context()
    {
        const auto& top_stack = context_stack.top();
        context_index.insert(top_stack.first, top_stack.second);
        context_stack.pop();
    }

public:
    /**
     * @brief The empty constructor
     */
    MutationEngine()
    {}

    /**
     * @brief A constructor
     * 
     * @param context_index is the genome context index
     * @param num_of_alleles is the number of alleles in wild-type cells
     * @param default_mutational_coefficients is the map of the default mutational signature coefficients
     * @param mutational_signatures is the map of the mutational signatures
     * @param mutational_properties is the species mutational properties
     * @param seed is the random generator seed
     */
    MutationEngine(ContextIndex<GENOME_WIDE_POSITION>& context_index,
                   const size_t& num_of_alleles,
                   const std::map<std::string, double>& default_mutational_coefficients,
                   const std::map<std::string, MutationalSignature>& mutational_signatures,
                   const SpeciesMutationalProperties& mutational_properties,
                   const int& seed=0):
        generator(seed), context_index(context_index), num_of_alleles(num_of_alleles),
        mutational_properties(mutational_properties)
    {
        for (const auto& [name, mutational_signature]: mutational_signatures) {
            inv_cumulative_SBSs[name] = get_inv_cumulative_SBS(mutational_signature);
        }

        add(0, default_mutational_coefficients);
    }

    /**
     * @brief Add a timed set of mutational signature coefficients
     * 
     * This method add a set mutational signature coefficients that will be applied from 
     * the specified simulation time on.
     * 
     * @param time is the simulation time from which the coefficients will be applied
     * @param mutational_coefficients is the map from the signature name to the coefficient
     * @return a reference to the updated object
     */
    MutationEngine& add(const Time& time, std::map<std::string, double>&& mutational_coefficients)
    {
        if (time<0) {
            throw std::domain_error("Simulation time is a non-negative value");
        }

        if (timed_mutational_coefficients.count(time)>0) {
            throw std::runtime_error("Another set of mutational coefficients has "
                                     "been set for the specified time");
        }

        for (const auto& [name, coeff]: mutational_coefficients) {
            if (inv_cumulative_SBSs.count(name)==0) {
                throw std::runtime_error("Unknown signature "+name);
            }
        }

        timed_mutational_coefficients.emplace(Time(time), std::move(mutational_coefficients));

        return *this;
    }

    /**
     * @brief Add a timed set of mutational signature coefficients
     * 
     * This method add a set mutational signature coefficients that will be applied from 
     * the specified simulation time on.
     * 
     * @param time is the simulation time from which the coefficients will be applied
     * @param mutational_coefficients is the map from the signature name to the coefficient
     * @return a reference to the updated object
     */
    MutationEngine& add(const Time& time, const std::map<std::string, double>& mutational_coefficients)
    {
        return add(time, std::map<std::string, double>(mutational_coefficients));
    }

    /**
     * @brief Place genomic mutations on a phylogenetic forest
     * 
     * @param forest is a phylogenetic forest
     * @param progress_bar is a progress bar pointer
     * @return the list of genomic mutations of the `forest`'s samples
     */
    std::list<SampleGenomeMutations> place_mutations(const Drivers::PhylogeneticForest& forest,
                                                     UI::ProgressBar *progress_bar=nullptr)
    {
        using namespace Races::Drivers::Simulation;
        using namespace Races::Passengers;

        GenomeMutations mutations(context_index, num_of_alleles);

        std::list<SampleGenomeMutations> sample_mutations;
        std::map<TissueSampleId, SampleGenomeMutations*> sample_mutation_map;
        for (const auto& sample: forest.get_samples()) {
            sample_mutations.push_back(SampleGenomeMutations(sample));
            sample_mutation_map.insert({sample.get_id(), &(sample_mutations.back())});
        }

        size_t visited_node = 0;
        for (const auto& root: forest.get_roots()) {
            place_mutations(sample_mutation_map, root, mutations, visited_node, progress_bar);
        }

        return sample_mutations;
    }

    /**
     * @brief Place genomic mutations on a phylogenetic forest
     * 
     * @param forest is a phylogenetic forest
     * @param progress_bar is a progress bar pointer
     * @return the list of genomic mutations of the `forest`'s samples
     */
    inline std::list<SampleGenomeMutations> place_mutations(const Drivers::PhylogeneticForest& forest,
                                                            UI::ProgressBar &progress_bar)
    {
        return place_mutations(forest, &progress_bar);
    }

    /**
     * @brief Bring back the context index to its original state
     * 
     * This method may be used after an exception to bring back the 
     * context index to its original state.
     */
    void reset_context_index()
    {
        while (!context_stack.empty()) {
            reinsert_last_context();
        }
    }
};

}   // Passengers

}   // Races

#endif // __RACES_MUTATION_ENGINE__
/**
 * @file read_simulator.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes to simulate sequencing
 * @version 0.5
 * @date 2023-10-26
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

#ifndef __RACES_READ_SIMULATOR__
#define __RACES_READ_SIMULATOR__

#include <list>
#include <map>
#include <string>
#include <sstream>
#include <filesystem>
#include <random>

#include "genomic_position.hpp"
#include "genome_mutations.hpp"

#include "fasta_reader.hpp"
#include "fasta_utils.hpp"

#include "progress_bar.hpp"
#include "variables.hpp"

#if WITH_MATPLOT
#include <matplot/matplot.h>
#endif // WITH_MATPLOT

namespace Races
{

namespace Passengers
{

/**
 * @brief The sequencing simulation namespace
 */

namespace SequencingSimulations
{

/**
 * @brief Statistics about simulated sequencing on a set of tissue samples
 */
struct Statistics
{
    /**
     * @brief Statistics about simulated sequencing on a tissue sample
     */
    struct SampleStatistics
    {
        std::map<ChromosomeId, std::vector<uint16_t>> coverage;  //!< The coverage profile
        std::map<SNV, size_t> SNV_occurrences;       //!< The SNV occurrences in the simulated reads

        /**
         * @brief Increase the coverage of a region
         * 
         * @param chromosome_id is the id of the chromosome containing the region whose coverage must be increase
         * @param begin_pos is the initial position of the region whose coverage must be increase
         * @param read_size is the size of the region whose coverage must be increase
         */
        void increase_coverage(const ChromosomeId& chromosome_id, const ChrPosition& begin_pos, const size_t& read_size);

        /**
         * @brief Increase the number of occurrences of an SNV
         * 
         * @param snv is the snv whose number of occurrences must be increased
         */
        void increase_SNV_occurrences(const SNV& snv);

        /**
         * @brief Add a chromosome to the statistics
         * 
         * @param chromosome_id is the id of the chromosome to be added
         * @param size is the size of the chromosome
         */
        inline void add_chromosome(const ChromosomeId& chromosome_id, const GenomicRegion::Length& size)
        {
            if (coverage.find(chromosome_id) == coverage.end()) {
                coverage.insert({chromosome_id, std::vector<uint16_t>(size,0)});
            }
        }

        /**
         * @brief Find the coverage in a genomic positions
         * 
         * @param position is the genomic position whose coverage is aimed
         * @return a constant reference to the coverage of the aimed position
         */
        const uint16_t& find_coverage(const GenomicPosition& position) const;
    };

private:
    /**
     * @brief Summarize the statistics of the sample
     * 
     * @return a sample statistics corresponding to the sums of the coverages
     *      and SNV occurrences of all the read simulated on the sample
     */
    SampleStatistics get_total() const;

public:

    std::map<std::string, SampleStatistics> sample_statistics;  //!< The statistics of a sample

    /**
     * @brief Get the statistics of a sample
     * 
     * @param sample_name is the aimed sample name
     * @return a reference to the sample statistics
     */
    inline SampleStatistics& operator[](const std::string& sample_name)
    {
        return sample_statistics[sample_name];
    }

    /**
     * @brief Save a csv file reporting the SNV VAFs
     * 
     * @param filename is the directory filename
     */
    void save_VAF_csv(const std::filesystem::path& filename) const;

#if WITH_MATPLOT
    /**
     * @brief Save a image representing the coverage of a chromosome
     * 
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
     *          coverage must be represented
     */
    void save_coverage_images(const std::filesystem::path& filename, const ChromosomeId& chromosome_id) const;

    /**
     * @brief Save a image representing the histogram of a chromosome SNVs
     * 
     * @param filename is the image filename
     * @param chromosome_id is the identifier of the chromosome whose
     *          SNVs must be represented
     */
    void save_SNV_histogram(const std::filesystem::path& filename, const ChromosomeId& chromosome_id) const;
#endif // WITH_MATPLOT
};

/**
 * @brief A read simulator 
 * 
 * @tparam RANDOM_GENERATOR is the random generator type used to sample genome
 */
template<typename RANDOM_GENERATOR=std::mt19937_64>
class ReadSimulator
{
public:
    /**
     * @brief The type of simulated reads
     */
    enum class ReadType {
        SINGLE_READ,
        PAIRED_READ
    };

    /**
     * @brief SAM generator mode
     */
    enum class Mode {
        CREATE,     // Create
        OVERWRITE,  // Overwrite files
        APPEND      // Append files
    };
private:
    RANDOM_GENERATOR random_generator;  //!< The random generator

    std::filesystem::path output_directory;         //!< The output directory
    std::filesystem::path ref_genome_filename;      //!< The reference genome FASTA filename

    ReadType read_type;     //!< The type of the produced reads

    size_t read_size;       //!< The produced-read size
    size_t insert_size;     //!< The paired-read insert size
    size_t num_of_reads;    //!< Number of already placed reads

    bool write_SAM;         //!< A Boolean flag to write SAM files

    /**
     * @brief A structure to maintain read simulation data
     */
    struct ReadSimulationData
    {
        size_t not_covered_allelic_size;   //!< Allelic size to cover yet
        size_t missing_templates;          //!< Number of templates to be placed

        /**
         * @brief The empty constructor
         */
        ReadSimulationData():
            ReadSimulationData(0,0)
        {}

        /**
         * @brief A constructor
         * 
         * @param not_covered_allelic_size is the allelic size to cover yet
         * @param missing_templates is the number of templates to be placed
         */
        ReadSimulationData(const size_t& not_covered_allelic_size,
                            const size_t& missing_templates):
            not_covered_allelic_size(not_covered_allelic_size), 
            missing_templates(missing_templates)
        {}
    };

    /**
     * @brief Chromosome data
     * 
     * This template represents chromosome data. The objects of this class store 
     * chromosome information, i.e., name, header, size, and chromosome identifier,
     * and, depending on the parameter, which must be in the hierarchy of the class
     * `Races::IO::FASTA::SequenceInfo`, may also maintain the chromosome nucleic
     * sequence.
     * 
     * @tparam DATA_TYPE is the base type of the template. It must be a inherited
     *      from `Races::IO::FASTA::SequenceInfo`
     */
    template<typename DATA_TYPE, std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    struct ChromosomeData : public DATA_TYPE
    {
        ChromosomeId chr_id;    //!< the chromosome id
    };

    struct FilterNonChromosomeSequence : public Races::IO::FASTA::SequenceFilter
    {
        ChromosomeId last_chr_id;   //!< the identifier of the last processed chromosome header

        inline bool operator()(const std::string& header)
        {
            return !Races::IO::FASTA::is_chromosome_header(header, last_chr_id);
        } 
    };

    /**
     * @brief Get the overall size of a list of sample mutations
     * 
     * This method computes the sum of the sizes of the genomic 
     * mutations in a list of sample mutations.
     * 
     * @param mutations_list is a list of sample mutations
     * @param progress_bar is a progress bar
     * @return the sum of the sizes of the genomic mutations in `mutations_list`
     */
    static size_t get_overall_size(const std::list<SampleGenomeMutations>& mutations_list,
                                   UI::ProgressBar& progress_bar)
    {
        progress_bar.set_message("Evaluating allelic size");
        size_t overall_size{0};
        for (const auto& sample_mutations: mutations_list) {
            for (const auto& mutations: sample_mutations.mutations) {
                overall_size += mutations.size();
                progress_bar.update_elapsed_time();
            }
        }

        return overall_size;
    }

    /**
     * @brief Get the overall allelic size for each element in a list of sample mutations
     * 
     * This method computes the sum of the allelic sizes of each sample genomic mutations 
     * in a list of sample mutations.
     * 
     * @param mutations_list is a list of sample mutations
     * @param progress_bar is a progress bar
     * @return a map associating the sample id to the sum of all the allelic sizes in 
     *      the sample mutations in `mutation_lists`
     */
    static std::map<Races::Drivers::Simulation::TissueSampleId, size_t>
    get_sample_allelic_sizes(const std::list<SampleGenomeMutations>& mutations_list,
                             UI::ProgressBar& progress_bar)
    {
        progress_bar.set_message("Evaluating allelic sizes");

        std::map<Races::Drivers::Simulation::TissueSampleId, size_t> allelic_sizes;
        for (const auto& sample_mutations: mutations_list) {
            auto& sample_allelic_size = allelic_sizes[sample_mutations.get_id()];
            sample_allelic_size = 0;
            for (const auto& mutations: sample_mutations.mutations) {
                sample_allelic_size += mutations.allelic_size();
                progress_bar.update_elapsed_time();
            }
        }

        return allelic_sizes;
    }

    /**
     * @brief Check whether all the genomic mutations in a list represent the same genome
     * 
     * This method checks whether all the genomic mutations in a list represent genomes having 
     * the same number of chromosomes and the chromosome sizes are consistent.
     * 
     * @tparam GENOME_MUTATION is the type of genome mutations
     * @param mutations is a list of genome mutations
     * @param progress_bar is a progress bar
     * @return `true` if and only if all the genomic mutations in `mutations` have the same 
     *          number of chromosomes and, for all pairs of elements in `mutations`, \f$m_1\f$
     *          and \f$m_2\f$, the $i$-th chromosomes in \f$m_1\f$ and \f$m_2\f$ are equal in
     *          size
     * @todo move this method in SampleGenomicMutations
     */
    template<typename GENOME_MUTATION, std::enable_if_t<std::is_base_of_v<Races::Passengers::GenomeMutations, GENOME_MUTATION>, bool> = true>
    static bool represents_same_genome(const std::list<GENOME_MUTATION>& mutations, UI::ProgressBar& progress_bar)
    {
        if (mutations.size()==0) {
            return true;
        }

        // get the first genome mutations in the list
        const auto& first_mutations = mutations.front();
        
        progress_bar.set_message("Checking input");
        // for all the genome mutations in the list
        for (const auto& cell_mutations : mutations) {

            // if the first element and the current one differ in the number of chromosome
            if (cell_mutations.get_chromosomes().size()!=first_mutations.get_chromosomes().size()) {
                return false;
            }

            // for all the chromosomes in the two genome mutations
            auto cell_it = cell_mutations.get_chromosomes().begin();
            auto first_it = first_mutations.get_chromosomes().begin();
            for (;cell_it != cell_mutations.get_chromosomes().end(); ++cell_it, ++first_it) {
                // if they differ because of their chromosome identifier or because of their 
                // sizes
                if ((cell_it->first != first_it->first)
                        ||(cell_it->second.size()!=first_it->second.size())) {
                    return false;
                }
            }
            progress_bar.update_elapsed_time();
        }

        // otherwise, return true
        return true;
    }

    /**
     * @brief Compute the CIGAR string of a genomic region potentially affected by SNVs
     * 
     * The CIGAR string a code that represents the alignment of a sequence over a reference
     * in SAM files (see [1]). This method considers a set of SNVs, a genomic position, 
     * and a size, and it produces the CIGAR code corresponding to a simulated read whose 
     * sequence correspond to that of the reference genome from the specified position with
     * the exception of positions in which the given SNVs were applied.
     * 
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification 
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
     * 
     * @param SNVs is a map from genomic positions to the corresponding SNVs
     * @param genomic_position is the initial position of the simulated read
     * @param read_size is the size of the simulated read
     * @return the CIGAR code corresponding to a simulated read having size `read_size` 
     *      whose sequence correspond to that of the sequence in the reference genome 
     *      from `genomic_position` with the exception of the genomic position in which 
     *      the SNVs in `SNVs` were placed
     */
    static std::string get_CIGAR(const std::map<GenomicPosition, SNV>& SNVs,
                                 const GenomicPosition& genomic_position, const size_t& read_size)
    {
        using namespace Races::Passengers;

        const ChrPosition last_position = genomic_position.position+read_size-1;

        std::ostringstream oss;

        // the first position in the last matching sequence
        size_t first_match = 0;

        // for all the SNVs placed on the simulated read
        for (auto snv_it = SNVs.lower_bound(genomic_position);
                snv_it != SNVs.end() && snv_it->second.position<=last_position; ++snv_it) {

            // get the SNV position in the simulated read
            const size_t SNV_pos = snv_it->second.position-genomic_position.position;

            // if we have matched some nucleotides
            if (SNV_pos+1>first_match) {
                // report the number of matched nucleotides
                oss << (SNV_pos-first_match) << "M";
            }

            // report the SNV mismatch
            oss << "1X";

            // the next matching sequence will start from the next 
            // position in the read
            first_match = SNV_pos+1;
        }

        // if the last position was a match
        if (genomic_position.position+first_match <= last_position+1) {
            // report the number of matched nucleotides
            oss << (last_position-genomic_position.position+1-first_match) << "M";
        }

        return oss.str();
    }

    /**
     * @brief Apply a set of SNVs to a fragment of the reference genome
     * 
     * @param nucleic_sequence is a fragment of the reference genome
     * @param SNVs is a map from genomic positions to the corresponding SNVs
     * @param genomic_position is the initial position of `nucleic_sequence` in
     *          the reference genome
     * @return a reference to the modified version of `nucleic_sequence` where 
     *          the SNVs that were placed in the reference genome region 
     *          corresponding to `nucleic_sequence` are applied
     */
    static std::string& apply_SNVs(std::string& nucleic_sequence,
                                   const std::map<GenomicPosition, SNV>& SNVs,
                                   const GenomicPosition& genomic_position)
    {
        using namespace Races::Passengers;

        const ChrPosition last_position = genomic_position.position+nucleic_sequence.size()-1;

        for (auto snv_it = SNVs.lower_bound(genomic_position);
                snv_it != SNVs.end() && snv_it->second.position<=last_position; ++snv_it) {
            const size_t SNV_pos = snv_it->second.position-genomic_position.position;

            {   // check the match between the SNV context and the corresponding nucleic 
                // triplet in `nucleic_sequence`

                // identify the first position in `nucleic_sequence` of the context
                size_t first_nucleotide = (SNV_pos==0?0:SNV_pos-1);

                // extract the `nucleic_sequence` nucleotides corresponding to the context
                std::string real_context = nucleic_sequence.substr(first_nucleotide, (SNV_pos==0?2:3));
                for (auto str_c: real_context) {
                    str_c = toupper(str_c);
                }

                // extract the part of the SNV context falling inside ``nucleic_sequence` 
                size_t length = real_context.size();
                std::string context = snv_it->second.context.get_sequence().substr((SNV_pos==0?1:0), length);

                // if they differ
                if (real_context!=context) {
                    std::ostringstream oss;

                    oss << "Context mismatch in chr. "<< GenomicPosition::chrtos(genomic_position.chr_id) 
                        << " in position " << snv_it->second.position << ": \"" 
                        << context << "\" expected and got \"" 
                        << real_context << "\"."<< std::endl 
                        << "Are you using a FASTA file different from that used for the context index?"
                        << std::endl;
                    throw std::domain_error(oss.str());
                }
            }

            // change the `nucleic_sequence` nucleotide according to the considered SNV
            nucleic_sequence[SNV_pos] = snv_it->second.mutated_base;
        }

        return nucleic_sequence;
    }

    /**
     * @brief Read the next chromosome data from the FASTA stream
     * 
     * @tparam DATA_TYPE if the type of the data to be read from the stream
     * @param FASTA_stream is the FASTA stream
     * @param chr_data is the object that will be filled with the chromosome data if
     *          some chromosome is read from `FASTA_stream`
     * @param progress_bar is the progress bar 
     * @return `true` if and only if a chromosome sequence is read from `FASTA_stream`
     */
    template<typename DATA_TYPE, std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    static bool read_next_chromosome_data(std::istream& FASTA_stream, ChromosomeData<DATA_TYPE>& chr_data, 
                                          UI::ProgressBar& progress_bar)
    {
        using namespace Races::IO::FASTA;

        FilterNonChromosomeSequence filter;

        if (DATA_TYPE::read(FASTA_stream, chr_data, filter, progress_bar)) {
            chr_data.chr_id = filter.last_chr_id;

            return true;
        }

        return false;
    }

    /**
     * @brief Get the template read data
     * 
     * @param template_first_position is the template first position
     * @param template_size is the template size
     * @return a vector containing a pair flag-first-base position for 
     *       each read in the template
     */
    std::vector<std::pair<int, ChrPosition>>
    get_template_read_data(const ChrPosition& template_first_position,
                           const size_t& template_size)
    {
        std::uniform_int_distribution<int> c_dist(0,1);

        bool rev_comp = c_dist(random_generator)==0; 

        std::vector<std::pair<int, ChrPosition>> template_read_data;

        if (read_type == ReadType::PAIRED_READ) {
            int flag = 0x1|0x2;
            template_read_data.push_back({flag|(rev_comp?0x10:0x0)|0x40,
                                          template_first_position});

            template_read_data.push_back({flag|(rev_comp?0x0:0x20)|0x80, 
                                          template_first_position+template_size-read_size});

            return template_read_data;
        }

        int flag = 0;
        if (rev_comp) {
            flag = 0x10;
        }
        
        return {{flag, template_first_position}};
    }

    /**
     * @brief Collect template statistics
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in] SNVs is a map from genomic positions to the corresponding SNVs
     * @param[in] read_first_position is the first position of the simulated read
     * @param[in] template_size is the size of the template
     * @param[in,out] sample_statistics are the statistics about a tissue sample
     */
    template<typename DATA_TYPE, 
         std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void collect_template_statistics(const ChromosomeData<DATA_TYPE>& chr_data,
                                     const std::map<GenomicPosition, SNV>& SNVs,
                                     const ChrPosition& template_first_position,
                                     const size_t& template_size, Statistics::SampleStatistics& sample_statistics)
    {
        const auto template_read_data = get_template_read_data(template_first_position, template_size);
        for (const auto& [flags, read_first_position] : template_read_data) {
        
            sample_statistics.increase_coverage(chr_data.chr_id, read_first_position, read_size);

            const GenomicPosition genomic_position{chr_data.chr_id, read_first_position};
            for (auto snv_it = SNVs.lower_bound(genomic_position);
                    snv_it != SNVs.end() && snv_it->second.position<read_first_position+read_size; ++snv_it) {
                sample_statistics.increase_SNV_occurrences(snv_it->second);
            }
        }
    }

    /**
     * @brief Write the SAM alignment line corresponding to a simulated read
     * 
     * This method writes in a stream the SAM alignment line corresponding to a simulated read
     * whose initial position on a reference genome chromosome, and SNVs are provided.
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param SAM_stream is the SAM file output stream
     * @param chr_data is the data about the chromosome from which the simulated read come from
     * @param SNVs is a map from genomic positions to the corresponding SNVs
     * @param read_first_position is the first position of the simulated read
     * @param template_size is the size of the template
     * @param sample_name is the tissue sample name
     */
    template<typename DATA_TYPE, 
         std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void write_SAM_template_alignments(std::ostream& SAM_stream, const ChromosomeData<DATA_TYPE>& chr_data,
                                       const std::map<GenomicPosition, SNV>& SNVs,
                                       const ChrPosition& template_first_position,
                                       const size_t& template_size, const std::string& sample_name="")
    {
        auto template_read_data = get_template_read_data(template_first_position, template_size);

        std::ostringstream oss_name;
        oss_name << "r" <<  std::setfill('0') << std::setw(11) 
                 << (num_of_reads/template_read_data.size());
        for (size_t i=0; i<template_read_data.size(); ++i) {
            const auto& read_first_position = template_read_data[i].second;

            const GenomicPosition genomic_position{chr_data.chr_id, read_first_position};

            std::string read;
            if constexpr(std::is_base_of_v<Races::IO::FASTA::Sequence, DATA_TYPE>) {
                read = chr_data.nucleotides.substr(read_first_position-1, read_size);
                apply_SNVs(read, SNVs, genomic_position);
            } else {
                read = "*";
            }

            size_t paired_idx = 1-i;

            std::string read_name = oss_name.str();
            
            if (read_type == ReadType::PAIRED_READ) {
                read_name += "/" + std::to_string(i+1);
            }

            const auto mapq = 33;

            SAM_stream << read_name                             // QNAME
                       << '\t' << template_read_data[i].first   // FLAG
                       << '\t' << chr_data.name                 // RNAME
                       << '\t' << read_first_position           // POS
                       << '\t' << mapq                          // MAPQ
                       << '\t' << get_CIGAR(SNVs, genomic_position, read_size); // CIGAR
            if (read_type == ReadType::PAIRED_READ) {
                SAM_stream << '\t' << oss_name.str() << "/" << (paired_idx+1)   // RNEXT
                           << '\t' << template_read_data[paired_idx].second     // PNEXT
                           << '\t' << (i==0?"":"-") << template_size;           // TLEN
            } else {
                SAM_stream << '\t' << '*'   // RNEXT
                           << '\t' << '0'   // PNEXT
                           << '\t' << '0';  // TLEN
            }
            SAM_stream << '\t' << read  // SEQ
                       << '\t' << '*';  // QUAL
            if (sample_name.size()!=0) {
                SAM_stream << '\t' << "RG:Z:" << sample_name;  // TAGS
            }
            SAM_stream << std::endl;
        }
    }

    /**
     * @brief Generate the reads in an allelic fragment
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param[in] fragment is the allele fragment for which reads must be generated 
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in,out] sample_simulation_data are the read simulation data relative to the considered `sample`
     * @param[in] sample_name is the name of the considered tissue sample
     * @param[in,out] sample_statistics are the statistics of the tissue sample identified by `sample_name`
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename DATA_TYPE=Races::IO::FASTA::Sequence,
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void generate_fragment_reads(const AlleleFragment& fragment, const ChromosomeData<DATA_TYPE>& chr_data,
                                 ReadSimulationData& sample_simulation_data, const std::string& sample_name, 
                                 Statistics::SampleStatistics& sample_statistics,
                                 const size_t& total_steps, size_t& steps,
                                 Races::UI::ProgressBar& progress_bar, std::ostream* SAM_stream)
    {
        auto template_size = ((read_type==ReadType::PAIRED_READ)
                              ?2*read_size+insert_size:read_size);

        if (fragment.size()>=template_size) {
            const double hit_probability = static_cast<double>(fragment.size())/
                                        sample_simulation_data.not_covered_allelic_size;
            std::binomial_distribution<size_t> b_dist(sample_simulation_data.missing_templates, 
                                                        hit_probability);

            size_t num_of_frag_templates = b_dist(random_generator);

            const double fragment_progress_ratio = (100*static_cast<double>(fragment.size())/num_of_frag_templates)/total_steps;

            const double current_progress = 100*static_cast<double>(steps)/total_steps;

            const auto& SNVs = fragment.get_SNVs();

            auto first_possible_begin = fragment.get_initial_position();
            auto last_possible_begin = static_cast<ChrPosition>(fragment.get_final_position())-read_size+1;

            std::uniform_int_distribution<ChrPosition> dist(first_possible_begin, last_possible_begin);
            size_t simulated_templates = 0;
            while (simulated_templates < num_of_frag_templates) {
                auto begin_pos = dist(random_generator);
                auto new_template_size = ((read_type==ReadType::PAIRED_READ)
                                            ?2*read_size+insert_size:read_size);
                auto end_pos = begin_pos + new_template_size - 1;
                if (end_pos > fragment.get_final_position() && begin_pos+1>new_template_size) {
                    end_pos = begin_pos;

                    begin_pos = end_pos+1-new_template_size-read_size; 
                }

                if (begin_pos >= fragment.get_initial_position()) {
                    if (SAM_stream != nullptr) {
                        write_SAM_template_alignments(*SAM_stream, chr_data, SNVs, begin_pos, 
                                                      new_template_size, sample_name);
                    }
                    collect_template_statistics(chr_data, SNVs, begin_pos,
                                                new_template_size, sample_statistics);
                    
                    num_of_reads += ((read_type == ReadType::PAIRED_READ)?2:1);
                    ++simulated_templates;
                    --sample_simulation_data.missing_templates;
                }
                progress_bar.set_progress(current_progress+simulated_templates*fragment_progress_ratio);
            }
        }

        sample_simulation_data.not_covered_allelic_size -= fragment.size();
        steps += fragment.size();
    }

    /**
     * @brief Generate simulated reads on a chromosome and write their SAM alignments
     * 
     * This method takes a list of genome mutations and a reference genome chromosome, it 
     * generates simulated reads over it up to a specified coverage of the mutated genomes, 
     * and write the corresponding SAM alignments in a stream.
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param[in] mutations_list is a list of sample mutations
     * @param[in] chr_data is the data about the chromosome from which the simulated read come from
     * @param[in,out] read_simulation_data are the read simulation data
     * @param[in] total_steps is the total number of steps required to complete the overall procedure
     * @param[in,out] steps is the number of performed steps
     * @param[in,out] progress_bar is the progress bar
     * @param[in,out] SAM_stream is the SAM file output stream
     */
    template<typename DATA_TYPE=Races::IO::FASTA::Sequence, 
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void generate_chromosome_reads(const std::list<SampleGenomeMutations>& mutations_list,
                                   const ChromosomeData<DATA_TYPE>& chr_data, 
                                   std::map<Races::Drivers::Simulation::TissueSampleId,
                                            ReadSimulationData>& read_simulation_data, 
                                   const size_t& total_steps, size_t& steps, 
                                   Races::UI::ProgressBar& progress_bar,
                                   std::ostream* SAM_stream=nullptr)
    {
        progress_bar.set_message("Processing chr. "+GenomicPosition::chrtos(chr_data.chr_id));

        Statistics statistics;

        for (const auto& sample_mutations: mutations_list) {
            const std::string& sample_name = sample_mutations.get_name();

            auto& sample_statistics = statistics[sample_name];
            auto& sample_simulation_data = read_simulation_data[sample_mutations.get_id()];

            for (const auto& cell_mutations: sample_mutations.mutations) {
                const auto& chr_mutations = cell_mutations.get_chromosome(chr_data.chr_id);

                sample_statistics.add_chromosome(chr_data.chr_id, chr_mutations.size());

                for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
                    for (const auto& [position, fragment] : allele.get_fragments()) {
                        generate_fragment_reads(fragment, chr_data, sample_simulation_data,
                                                sample_name, sample_statistics, total_steps,
                                                steps, progress_bar, SAM_stream);
                    }

                    progress_bar.set_progress(100*steps/total_steps);
                }
            }
        }

        progress_bar.set_message("Saving statistics");

        const std::string base_filename = "chr_"+GenomicPosition::chrtos(chr_data.chr_id);
        statistics.save_VAF_csv(output_directory/(base_filename+"_VAF.csv"));
#if WITH_MATPLOT
        progress_bar.update_elapsed_time();

        progress_bar.set_message("Saving "+base_filename+" images");
        statistics.save_coverage_images(output_directory/(base_filename+"_coverage.jpg"), chr_data.chr_id);
        progress_bar.update_elapsed_time();
        statistics.save_SNV_histogram(output_directory/(base_filename+"_hist.jpg"), chr_data.chr_id);
#endif // WITH_MATPLOT
    }

    /**
     * @brief Get the SAM stream
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param chr_data in the information available about the chromosome to be considered
     * @param mutations_list is a list of sample mutations
     * @return the SAM stream
     */
    template<typename DATA_TYPE=Races::IO::FASTA::Sequence,
            std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    std::ofstream
    get_SAM_stream(const ChromosomeData<DATA_TYPE>& chr_data,
                   const std::list<SampleGenomeMutations>& mutations_list) const
    {
        auto chr_str = GenomicPosition::chrtos(chr_data.chr_id);

        std::string SAM_filename = "chr_"+chr_str+".sam";
        auto SAM_file_path = output_directory/SAM_filename;

        if (std::filesystem::exists(SAM_file_path)) {
            if (!std::filesystem::is_regular_file(SAM_file_path)) {
                throw std::runtime_error("\""+std::string(SAM_file_path)+"\" is not a regular file");
            }

            return std::ofstream(SAM_file_path, std::ofstream::app);
        }
    
        auto SAM_stream = std::ofstream(SAM_file_path);

        // write the SAM header
        SAM_stream << "@HD\tVN:1.6\tSO:unknown" << std::endl;
        SAM_stream << "@SQ\tSN:" << chr_data.name << "\tLN:" << chr_data.length 
                << "\tAN:chromosome" << chr_str << ",chr" << chr_str 
                << ",chromosome_" << chr_str << ",chr_" << chr_str << std::endl;
        for (const auto& sample_mutations : mutations_list) {
            SAM_stream << "@RG\tID:" << sample_mutations.get_name() << std::endl;
        }

        return SAM_stream;
    }

    /**
     * @brief Compute the initial read simulation data
     * 
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @return the initial read simulation data for each sample
     */
    std::map<Races::Drivers::Simulation::TissueSampleId, ReadSimulationData>
    get_initial_read_simulation_data(const std::list<SampleGenomeMutations>& mutations_list,
                                     const double& coverage) const
    {
        std::map<Races::Drivers::Simulation::TissueSampleId, ReadSimulationData> read_simulation_data;

        size_t total_read_size = (read_type==ReadType::PAIRED_READ?2:1)*read_size;
        for (const auto& sample_mutations : mutations_list) {
            auto& sample_data = read_simulation_data[sample_mutations.get_id()];
            if (sample_mutations.mutations.size()>0) {
                sample_data.missing_templates = 
                                 static_cast<size_t>((sample_mutations.mutations.front().size()*
                                                     coverage)/total_read_size);
            }
            for (const auto& cell_mutations: sample_mutations.mutations) {
                sample_data.not_covered_allelic_size += cell_mutations.allelic_size();
            }
        }

        return read_simulation_data;
    }

    /**
     * @brief Get the number of chromosomes in a genome
     * 
     * @param mutations_list is a list of sample mutations
     * @return the number of chromosomes in a genome 
     */
    static size_t get_number_of_chromosomes_in_a_genome(const std::list<SampleGenomeMutations>& mutations_list)
    {
        for (const auto& sample_mutations : mutations_list) {
            for (const auto& mutations : sample_mutations.mutations) {
                return mutations.get_chromosomes().size();
            }
        }

        return 0;
    }

    /**
     * @brief Generate simulated reads
     * 
     * This method takes a list of sample mutations and it generates simulated reads up to a 
     * specified coverage of the mutated genomes. The base-coverage and SNV occurrences are collected
     * and, upon request, the SAM alignments corresponding to the produced reads are saved. 
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @param progress_bar is the progress bar
     */
    template<typename DATA_TYPE=Races::IO::FASTA::Sequence, 
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void generate_reads(const std::list<SampleGenomeMutations>& mutations_list,
                        const double& coverage, UI::ProgressBar& progress_bar)
    {
        std::ifstream ref_stream(ref_genome_filename);

        auto read_simulation_data = get_initial_read_simulation_data(mutations_list, coverage);

        size_t steps{0};
        const auto num_of_chr_in_a_genome = get_number_of_chromosomes_in_a_genome(mutations_list);
        size_t total_steps = 2*num_of_chr_in_a_genome*mutations_list.size();
        for (const auto& [sample_id, sample_data]: read_simulation_data) {
            total_steps += sample_data.not_covered_allelic_size;
        }

        ChromosomeData<DATA_TYPE> chr_data;
        
        progress_bar.set_progress((100*steps)/total_steps, "Reading next chromosome");
        while (read_next_chromosome_data(ref_stream, chr_data, progress_bar)) {
            progress_bar.set_progress((100*(++steps))/total_steps);

            if (write_SAM) {
                std::ofstream SAM_stream = get_SAM_stream(chr_data, mutations_list);

                generate_chromosome_reads(mutations_list, chr_data, read_simulation_data, 
                                          total_steps, steps, progress_bar, &SAM_stream);
            } else {
                generate_chromosome_reads(mutations_list, chr_data, read_simulation_data,
                                          total_steps, steps, progress_bar);
            }
        }
        
        progress_bar.set_message("Read simulated");
        progress_bar.set_progress(100);
    }

    /**
     * @brief A constructor
     * 
     * @param output_directory is the output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_type is the type of the produced-read, i.e., single or paired-read
     * @param read_size is the size of the output reads
     * @param insert_size is the size of the insert
     * @param mode is the SAM generator output mode
     * @param seed is the random generator seed
     */
    ReadSimulator(const std::string& output_directory, const std::string& ref_genome_filename,
                  const ReadType read_type, const size_t& read_size, const size_t& insert_size, 
                  const Mode mode, const int& seed):
        random_generator(seed), output_directory(output_directory), ref_genome_filename(ref_genome_filename), 
        read_type(read_type), read_size(read_size), insert_size(insert_size), num_of_reads(0),
        write_SAM(false)
    {
        namespace fs = std::filesystem;

        if (fs::exists(this->output_directory)) {
            if (mode==Mode::CREATE) {
                throw std::domain_error("\""+output_directory+"\" already exists");
            }

            if (!fs::is_directory(this->output_directory)) {
                throw std::domain_error("\""+output_directory+"\" is not a directory");
            }

            if (mode==Mode::OVERWRITE) {
                fs::remove_all(this->output_directory);
            }
        }

        fs::create_directory(this->output_directory);

        if (!fs::exists(this->ref_genome_filename)) {
            throw std::domain_error("\""+ref_genome_filename+"\" does not exist");
        }

        if (!fs::is_regular_file(this->ref_genome_filename)) {
            throw std::domain_error("\""+ref_genome_filename+"\" is not a regular file");
        }

        std::ifstream ref_stream(ref_genome_filename);
        char c = ref_stream.get();

        if (c != EOF && c != '>') {
            throw std::domain_error("\""+ref_genome_filename+"\" is not a FASTA file");
        }
    }
public:

    /**
     * @brief The empty constructor
     */
    ReadSimulator():
        random_generator(), output_directory(""), ref_genome_filename(""), 
        read_type(ReadType::SINGLE_READ), read_size(0), insert_size(0), num_of_reads(0),
        write_SAM(false)
    {}

    /**
     * @brief Create a simulator that produces single reads
     * 
     * @param output_directory is the output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_size is the size of the output reads
     * @param mode is the SAM generator output mode
     * @param seed is the random generator seed
     */
    ReadSimulator(const std::string& output_directory, const std::string& ref_genome_filename,
                  const size_t& read_size, const Mode mode=Mode::CREATE, const int& seed=0):
        ReadSimulator(output_directory, ref_genome_filename, ReadType::SINGLE_READ, read_size,
                      0, mode, seed)
    {}

    /**
     * @brief Create a simulator that produces paired reads
     * 
     * @param output_directory is the SAM output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_size is the size of the output reads
     * @param insert_size is the size of the insert
     * @param mode is the SAM generator output mode
     * @param seed is the random generator seed
     */
    ReadSimulator(const std::string& output_directory, const std::string& ref_genome_filename,
                  const size_t& read_size, const size_t& insert_size, const Mode mode=Mode::CREATE, 
                  const int& seed=0):
        ReadSimulator(output_directory, ref_genome_filename, ReadType::PAIRED_READ, read_size,
                      insert_size, mode, seed)
    {}

    /**
     * @brief Generate simulated reads for a list of genome mutations
     * 
     * This method generates simulated reads for a set of cell mutated genomes and, 
     * if requested, writes the corresponding SAM alignments. The number of simulated reads 
     * depends on the specified coverage.
     * 
     * @param mutations_list is a list of sample mutations
     * @param coverage is the aimed coverage
     * @param save_SAM is a Boolean flag to save the simulated read in a file
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     * @return a reference to the updated object
     */
    ReadSimulator& operator()(const std::list<Races::Passengers::SampleGenomeMutations>& mutations_list,
                              const double& coverage, const bool with_sequences=true,
                              const bool quiet=false)
    {
        using namespace Races;

        if (coverage < 0) {
            std::ostringstream oss;

            oss << "coverage must be non-negative. Got " << coverage;
            throw std::domain_error(oss.str());
        }

        if (mutations_list.size()==0) {
            return *this;
        }

        UI::ProgressBar progress_bar(quiet);

        if (with_sequences) {
            generate_reads<Races::IO::FASTA::Sequence>(mutations_list, coverage, progress_bar);
        } else {
            generate_reads<Races::IO::FASTA::SequenceInfo>(mutations_list, coverage, progress_bar);
        }

        return *this;
    }

    /**
     * @brief Enable SAM file writing
     * 
     * @param enable is a flag to enable(true)/disable(false) SAM file writing
     */
    inline void enable_SAM_writing(const bool& enable)
    {
        write_SAM = enable;
    }
};

}   // SequencingSimulations

}   // Passengers

}   // Races

#endif // __RACES_READ_SIMULATOR__
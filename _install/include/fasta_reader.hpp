/**
 * @file fasta_reader.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a FASTA file reader and support structures
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

#ifndef __RACES_FASTA_READER__
#define __RACES_FASTA_READER__

#include <string>
#include <istream>
#include <algorithm>

#include "progress_bar.hpp"

namespace Races
{

namespace IO
{

/**
 * @brief A name space for FASTA file format
 */
namespace FASTA
{

/**
 * @brief A class to filter FASTA sequences
 * 
 * The objects of this class filter FASTA sequences according to 
 * their sequence headers. The class offers the method 
 * `operator()(const std::string& header)` that returns `true` 
 * when the sequence must be filtered and `false` otherwise.
 */
struct SequenceFilter
{
    /**
     * @brief The filtering method
     * 
     * This method returns `true` when the sequence must be 
     * filtered and `false` otherwise according to the sequence
     * header.
     * 
     * @param header is the FASTA sequence header
     * @return `true` if and only if the sequence must be filtered
     */
    inline constexpr bool operator()(const std::string& header) const
    {
        (void)header;
        
        return false;
    } 
};

/**
 * @brief A class to represent FASTA sequence information
 */
class SequenceInfo
{
    /**
     * @brief Read the next content of the stream up-to a new sequence
     * 
     * @param fasta_stream is an input stream referring to a FASTA file
     * @param progress_bar is a progress bar
     */
    static void filter_remaining_sequence(std::istream& fasta_stream, UI::ProgressBar& progress_bar);

protected:
    /**
     * @brief Read FASTA sequence information from a stream
     * 
     * @tparam FILTER is the type of sequence filter
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] seq_info is the object that will be filled by the read information
     * @param[out] nucleotides is a pointer to the string that will be filled by the sequence nucleotides
     * @param[in] filter is the filter selecting the appropriate FASTA sequences
     *              according to their headers
     * @param[in,out] progress_bar is a progress bar
     * @return `true` if and only if a sequence that is not filtered by `filter`
     *           has been read from `fasta_stream`. If the method returns `true`, 
     *           then `seq_info` and, whenever `nucleotides` is not `nullptr`, 
     *           `*nucleotides` are updated according to the read values
     */
    template<typename FILTER, std::enable_if_t<std::is_base_of_v<SequenceFilter, FILTER>,bool> = true>
    static bool read(std::istream& fasta_stream, SequenceInfo& seq_info, std::string* nucleotides,
                     FILTER& filter, UI::ProgressBar& progress_bar)
    {
        char c = fasta_stream.get();
        while (c != EOF) {
            if (c!='>') {
                fasta_stream.unget();
                std::string char_s(">");

                char_s[0] = c;
                throw std::runtime_error("expecting '>' got '"+char_s+"'");
            }
            std::string header;
            getline(fasta_stream, header);

            auto seq_name = header.substr(0, header.find(" "));

            progress_bar.set_message("Found "+seq_name);
            if (!filter(header)) {
                progress_bar.set_message("Reading "+seq_name);

                seq_info.length=0;
                seq_info.header = std::move(header);
                seq_info.name = std::move(seq_name);
                if (nucleotides != nullptr) {
                    *nucleotides = "";
                }

                size_t counter{0};
                while ((c = fasta_stream.get()) != EOF) {
                    std::string line;
                    fasta_stream.unget();
                    if (c=='>') {
                        return true;
                    }

                    getline(fasta_stream, line);

                    // let time elapse in progress bar
                    if ((counter = (counter+1)%1000) == 0) {
                        progress_bar.set_progress(progress_bar.get_progress());
                    }

                    line.erase(remove(line.begin(), line.end(), ' '), line.end());
                    line.erase(remove(line.begin(), line.end(), '\t'), line.end());
                    seq_info.length += line.size();
                    if (nucleotides != nullptr) {
                        nucleotides->append(line);
                    }
                }
                return true;
            } else {
                filter_remaining_sequence(fasta_stream, progress_bar);

                c = fasta_stream.get();
            }
        }

        return false;
    }

public:
    std::string name;   //!< The name of the FASTA sequence
    std::string header; //!< The header of the FASTA sequence
    size_t length;      //!< The length of the FASTA sequence

    /**
     * @brief Read FASTA sequence information from a stream
     * 
     * @param fasta_stream is a stream referring to a FASTA file
     * @param seq_info is the object that will be filled by the read information
     * @param progress_bar is a progress bar
     * @return `true` if and only if a sequence has been read from `fasta_stream`. 
     *           If the method returns `true`, then `seq_info` is updated
     *           according to the read values
     */
    static bool read(std::istream& fasta_stream, SequenceInfo& seq_info, UI::ProgressBar& progress_bar);

    /**
     * @brief Read FASTA sequence information from a stream
     * 
     * @tparam FILTER is the type of sequence filter
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] seq_info is the object that will be filled by the read information
     * @param[in] filter is the filter selecting the appropriate FASTA sequences
     *              according to their headers
     * @param[in,out] progress_bar is a progress bar
     * @return `true` if and only if a sequence that is not filtered by `filter`
     *           has been read from `fasta_stream`. If the method returns `true`, 
     *           then `seq_info` is updated according to the read values
     */
    template<typename FILTER, std::enable_if_t<std::is_base_of_v<SequenceFilter, FILTER>,bool> = true>
    inline static bool read(std::istream& fasta_stream, SequenceInfo& seq_info, FILTER& filter,
                            UI::ProgressBar& progress_bar)
    {
        return read(fasta_stream, seq_info, nullptr, filter, progress_bar);
    }

    /**
     * @brief Read FASTA sequence information from a stream
     * 
     * @tparam FILTER is the type of sequence filter
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] seq_info is the object that will be filled by the read information
     * @return `true` if and only if a sequence has been read from `fasta_stream`.
     *           If the method returns `true`, then `seq_info` is updated 
     *           according to the read values
     */
    static bool read(std::istream& fasta_stream, SequenceInfo& seq_info);

    /**
     * @brief Read FASTA sequence information from a stream
     * 
     * @tparam FILTER is the type of sequence filter
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] seq_info is the object that will be filled by the read information
     * @param[in] filter is the filter selecting the appropriate FASTA sequences
     *              according to their headers
     * @return `true` if and only if a sequence that is not filtered by `filter`
     *           has been read from `fasta_stream`. If the method returns `true`, 
     *           then `seq_info` is updated according to the read values
     */
    template<typename FILTER, std::enable_if_t<std::is_base_of_v<SequenceFilter, FILTER>,bool> = true>
    inline static bool read(std::istream& fasta_stream, SequenceInfo& seq_info, FILTER& filter)
    {
        UI::ProgressBar progress_bar(true);

        return SequenceInfo::read(fasta_stream, seq_info, filter, progress_bar);
    }
};

/**
 * @brief A class to represent FASTA sequence
 */
struct Sequence : public SequenceInfo
{
    std::string nucleotides;    //!< The string of sequence nucleotides

    /**
     * @brief Read FASTA sequence from a stream
     * 
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] sequence is the object that will be filled by the read sequence
     * @param[in,out] progress_bar is a progress bar
     * @return `true` if and only if a sequence has been read from `fasta_stream`. 
     *           If the method returns `true`, then `sequence` is updated
     *           according to the read values
     */
    static bool read(std::istream& fasta_stream, Sequence& sequence, UI::ProgressBar& progress_bar);

    /**
     * @brief Read FASTA sequence from a stream
     * 
     * @tparam FILTER is the type of sequence filter
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] sequence is the object that will be filled by the read sequence
     * @param[in] filter is the filter selecting the appropriate FASTA sequences
     *              according to their headers
     * @param[in,out] progress_bar is a progress bar
     * @return `true` if and only if a sequence that is not filtered by `filter`
     *           has been read from `fasta_stream`. If the method returns `true`, 
     *           then `seq_info` is updated according to the read values
     */
    template<typename FILTER, std::enable_if_t<std::is_base_of_v<SequenceFilter, FILTER>,bool> = true>
    inline static bool read(std::istream& fasta_stream, Sequence& sequence, FILTER& filter,
                            UI::ProgressBar& progress_bar)
    {
        return SequenceInfo::read(fasta_stream, sequence, &(sequence.nucleotides), filter, progress_bar);
    }

    /**
     * @brief Read FASTA sequence from a stream
     * 
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] sequence is the object that will be filled by the read sequence
     * @return `true` if and only if a sequence has been read from `fasta_stream`. 
     *           If the method returns `true`, then `sequence` is updated
     *           according to the read values
     */
    static bool read(std::istream& fasta_stream, Sequence& sequence);

    /**
     * @brief Read FASTA sequence from a stream
     * 
     * @tparam FILTER is the type of sequence filter
     * @param[in,out] fasta_stream is a stream referring to a FASTA file
     * @param[out] sequence is the object that will be filled by the read sequence
     * @param[in] filter is the filter selecting the appropriate FASTA sequences
     *              according to their headers
     * @return `true` if and only if a sequence that is not filtered by `filter`
     *           has been read from `fasta_stream`. If the method returns `true`, 
     *           then `seq_info` is updated according to the read values
     */
    template<typename FILTER, std::enable_if_t<std::is_base_of_v<SequenceFilter, FILTER>,bool> = true>
    inline static bool read(std::istream& fasta_stream, Sequence& sequence, FILTER& filter)
    {
        UI::ProgressBar progress_bar(true);

        return Sequence::read(fasta_stream, sequence, filter, progress_bar);
    }
};

}   // FASTA

}   // IO

}   // Races

#endif // __RACES_FASTA_READER__
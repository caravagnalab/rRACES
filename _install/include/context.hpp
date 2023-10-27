/**
 * @file context.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines mutational contexts and extended context automata
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

#ifndef __RACES_CONTEXT__
#define __RACES_CONTEXT__

#include <cstdint>
#include <functional>   // std::less
#include <string>
#include <iostream>
#include <array>

#include "archive.hpp"

namespace Races 
{

namespace Passengers
{

/**
 * @brief A class to represent mutational context
 * 
 * A context is a nucleic triplet where a mutation may 
 * occur.
 * Every context univocally corresponds to a code that 
 * represents the triplet. This class guarantees the 
 * mutational context code to be in the interval 
 * natural [0,63]: 5', 3', and the central nucleotide 
 * have 4 possible values, i.e., A, C, G, and T.
 * Moreover, all the codes for contexts having either 
 * a 'C' and a 'T' as central nucleotide are guaranteed
 * to be in the interval [0,31].
 */
struct MutationalContext
{
    /**
     * @brief The type of mutational context type
     */
    using CodeType = uint8_t;

private:
    CodeType code;   //!< the code associated to the context
public:
    /**
     * @brief The empty constructor
     */
    MutationalContext();
    
    /**
     * @brief A constructor
     * 
     * @param nucleic_triplet is the string of the nucleic triplet
     */
    MutationalContext(const char* nucleic_triplet);

    /**
     * @brief A constructor
     * 
     * @param nucleic_triplet is the string of the nucleic triplet
     */
    MutationalContext(const std::string& nucleic_triplet);

    /**
     * @brief Get the mutational context code 
     * 
     * @return a constant reference to the mutational context code 
     */
    inline const CodeType& get_code() const
    {
        return code;
    }

    /**
     * @brief Get the nucleic sequence of the context
     * 
     * @return the nucleic sequence of the context
     */
    std::string get_sequence() const;

    /**
     * @brief Get the context central nucleotide
     * 
     * @return the context central nucleotide
     */
    char get_central_nucleotide() const;

    /**
     * @brief Get the complement mutational context
     * 
     * @return the complement mutational context 
     */
    MutationalContext get_complement() const;

    /**
     * @brief Test whether two mutational contexts are equivalent
     * 
     * @param context is the mutational context to compare
     * @return `true` if and only if the two mutational contexts represent
     *      the same nucleic triplet
     */
    inline bool operator==(const MutationalContext& context) const
    {
        return get_code() == context.get_code();
    }

    /**
     * @brief Test whether two mutational contexts differ
     * 
     * @param context is the mutational context to compare
     * @return `true` if and only if the two mutational contexts represent
     *      different nucleic triplets
     */
    inline bool operator!=(const MutationalContext& context) const
    {
        return get_code() != context.get_code();
    }

    /**
     * @brief Save a mutational context in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & code;
    }

    /**
     * @brief Load a mutational context from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded mutational context
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static MutationalContext load(ARCHIVE& archive)
    {
        MutationalContext mcode;

        archive & mcode.code;

        return mcode;
    }

    /**
     * @brief Get the code of the complement mutational context
     * 
     * @param code is the mutational context code whose complement is request 
     * @return the code of the complement context of the 
     *      context whose code is `code` 
     */
    static CodeType get_complement(const CodeType& code);

    /**
     * @brief Get the complement base
     * 
     * @param base is the base whose complement is request 
     * @return the complement base of `base` 
     */
    static char get_complement(const char& base);

    /**
     * @brief Get the complement sequence
     * 
     * @param sequence is the sequence whose complement is request 
     * @return the complement sequence of `sequence` 
     */
    static std::string get_complement(const std::string& sequence);

    /**
     * @brief Check whether a symbol represents a nucleic base
     * 
     * @param symbol is the symbol to be tested 
     * @return `true` if and only if symbol represents a base, i.e., 
     *       it is one among 'A', 'C', 'G', 'T', 'a', 'c', 'g', or 't'.
     */
    static bool is_a_base(const char& symbol);
};


/**
 * @brief An automaton for parsing extended context
 * 
 * An extended context is a triplet in the alphabet 
 * 'A', 'C', 'G', 'T', and 'N'.
 * An object of this class is an automaton whose state 
 * correspond to an extended context. The automaton state 
 * always matches the extended context of the last read
 * nucleotides assuming that "NNN" is the initial
 * configuration. 
 */
struct ExtendedContextAutomaton
{
    /**
     * @brief The type of code for bases
     */
    using BaseCodeType = uint8_t;
    
    /**
     * @brief The type of the edges from a node
     * 
     * This type represents edges leaving a node: the key in the 
     * map is the code of the newly read base; the value is the 
     * new context code.
     */
    using FromNodeEdgeType = std::map<BaseCodeType, MutationalContext::CodeType>;

    std::array<bool, 125> is_a_context;             //!< a Boolean vector to flag proper context 
    std::array<MutationalContext, 125> contexts;    //!< the state corresponding contexts
    std::array<FromNodeEdgeType, 125> edges;        //!< the automata edges

    MutationalContext::CodeType state;  //!< the current automaton state

    /**
     * @brief Get a code for a character
     * 
     * @param character is a possible base
     * @return a number in the interval [0,5]. The 
     *      value 5 is returned for all the characters 
     *      different from  'A', 'C', 'G', 'T', and 'N'.
     */
    static BaseCodeType base2code(const char& character);

    /**
     * @brief Get the state corresponding to an extended context
     * 
     * @param first is the first character of the extended context
     * @param second is the second character of the extended context
     * @param third is the third character of the extended context
     * @return the automaton state corresponding to the extended
     *      context formed by the three characters.
     */
    static MutationalContext::CodeType get_state_for(const char& first, const char& second,
                                                     const char& third);

public:
    /**
     * @brief Build an extended context automaton
     */
    ExtendedContextAutomaton();

    /**
     * @brief Get the current automaton state
     * 
     * @return return the current automaton state
     */
    inline const MutationalContext::CodeType& get_state() const
    {
        return state;
    }

    /**
     * @brief Check whether the automaton state corresponds to a proper context
     * 
     * @return `true` if and only if the last three nucleotides read by the 
     *      automaton are among 'A', 'C', 'G', 'T'
     */
    inline const bool& read_a_context() const
    {
        return is_a_context[state];
    }

    /**
     * @brief Get the context corresponding to the current state
     * 
     * @return the context corresponding to the current state
     */
    inline const MutationalContext& get_context() const
    {
        return contexts[state];
    }

    /**
     * @brief Read a character and update the automaton state
     * 
     * @param character is the read character
     * @return `true` if and only if the read character is a valid 
     *      nucleotide character, i.e., it is among 'A', 'C', 'G',
     *      'T', and 'N'.
     */
    bool update_state(const char& character);

    /**
     * @brief Reset the automaton state
     * 
     * @return a reference to the updated automaton
     */
    ExtendedContextAutomaton& reset();

    /**
     * @brief Check whether the character is a nucleotide
     * 
     * @param character is a character
     * @return `true` if and only if the read character is a valid 
     *      nucleotide character, i.e., it is among 'A', 'C', 'G',
     *      'T', and 'N'.
     */
    static inline bool is_a_nucleotide(const char& character)
    {
        return base2code(character)<5;
    }
};

}   // Passengers

}   // Races

template<>
struct std::less<Races::Passengers::MutationalContext>
{
    inline bool operator()(const Races::Passengers::MutationalContext &lhs,
                           const Races::Passengers::MutationalContext &rhs) const
    {
        return lhs.get_code() < rhs.get_code();
    }
};

namespace std
{

/**
 * @brief Stream the mutational context in a stream
 * 
 * @param out is the output stream
 * @param context is the mutational context to stream
 * @return a reference to the output stream
 */
inline std::ostream& operator<<(std::ostream& out, const Races::Passengers::MutationalContext& context)
{
    return (out << context.get_sequence());
}

}   // std

#endif // __RACES_CONTEXT__
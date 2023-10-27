/**
 * @file snv_signature.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines Single Variation Mutation mutational signature
 * @version 0.12
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

#ifndef __RACES_SNV_DISTRIBUTION__
#define __RACES_SNV_DISTRIBUTION__

#include <string>
#include <map>
#include <set>
#include <functional> // std::less
#include <iostream>
#include <sstream>

#include "context.hpp"

namespace Races 
{

namespace Passengers
{

/**
 * @brief A class to represent mutational type
 * 
 * A mutational type is a mutational context and the 
 * replacing nucleic base for the central nucleotide.
 * For consistency with the literature, the mutational 
 * context will always be either a 'C' or a 'T'. We 
 * call it normal context. If the provided context 
 * is not a normal context, the complement context 
 * and replacing nucleic base are used.
 */
class MutationalType
{
    MutationalContext context;    //!< the normal context
    char replace_base;            //!< the replace base
public:
    /**
     * @brief The empty constructor
     */
    MutationalType();

    /**
     * @brief A constructor
     * 
     * @param context is the mutational type context
     * @param replace_base is the base replacing the context central nucleotide
     */
    MutationalType(const MutationalContext& context, const char& replace_base);

    /**
     * @brief A constructor
     * 
     * @param context is the mutational type context
     * @param replace_base is the base replacing the context central nucleotide
     */
    MutationalType(const std::string& context, const char& replace_base);

    /**
     * @brief A constructor
     * 
     * A mutational type is conventionally represented by a string in the 
     * form `X[Y>W]K` where `X` and `K` and the bases on 5' and 3', respectively,  
     * and `Y` and `W` are the central nucleotide before and after the mutation.
     * This constructor takes as parameter a string in this format and build 
     * the corresponding `MutationalType` object.
     * 
     * @param type is the textual representation of a mutational type
     */
    explicit MutationalType(const std::string& type);

    /**
     * @brief Get the mutational context
     * 
     * @return the mutational context
     */
    inline const MutationalContext& get_context() const
    {
        return context;
    }

    /**
     * @brief Get the replace base
     * 
     * @return the replace base
     */
    inline const char& get_replace_base() const
    {
        return replace_base;
    }

    /**
     * @brief Get complement base of the replace base
     * 
     * @return the replace base
     */
    inline char get_complement_replace_base() const
    {
        return MutationalContext::get_complement(replace_base);
    }

    /**
     * @brief Test whether two mutational types are equivalent
     * 
     * @param type is the mutational type to compare
     * @return `true` if and only if the two mutational types refer
     *      to equivalent contexts and the same replace base
     */
    inline bool operator==(const MutationalType& type) const
    {
        return (context == type.context) && (replace_base == type.replace_base);
    }

    /**
     * @brief Test whether two mutational types differ
     * 
     * @param type is the mutational type to compare
     * @return `true` if and only if the two mutational types refer
     *      to different contexts or different replace bases
     */
    inline bool operator!=(const MutationalType& type) const
    {
        return !(*this == type);
    }
};

}   // Passengers

}   // Races

template<>
struct std::less<Races::Passengers::MutationalType>
{
    bool operator()(const Races::Passengers::MutationalType &lhs,
                    const Races::Passengers::MutationalType &rhs) const;
};

namespace Races 
{

namespace Passengers
{


class MutationalSignature;

/**
 * @brief A class to represent the result of a mutational signature expression
 * 
 * This class is meant to represent temporary object that are evaluated during 
 * the computation of expressions of the kind:
 *  
 * \f$\alpha_1 * \beta_1 + alpha_2 * beta_2 + \ldots\f$
 * 
 * where the \f$\alpha_i\f$'s are real values in the interval \f$[0,1]\f$ and 
 * the \f$\beta_i\f$'s are `MutationalSignature` objects. 
 * 
 * Even if the final result of above expression is a mutational signature, the
 * partial results may be different from a probability distribution and that 
 * is why this class is needed.
 */
class MutationalSignatureExprResult
{
    std::map<MutationalType, double> value_map; //!< the mutational type-value map

    /**
     * @brief The constructor
     * 
     * This constructor is private and it is meant to be exclusively called by 
     * `MutationalSignature`'s methods.
     * 
     * @param value_map is a mutational type-value map
     */
    MutationalSignatureExprResult(const std::map<MutationalType, double>& value_map);
public:
    /**
     * @brief The empty constructor
     */
    MutationalSignatureExprResult();

    /**
     * @brief Cast to `MutationalSignature`
     * 
     * This method tries to cast a mutational signature expression to a 
     * mutational signature. When the expression does not represent a 
     * probability distribution a `std::domain_error` is thrown.
     * 
     * @return the corresponding `MutationalSignature` object
     */
    operator MutationalSignature();

    /**
     * @brief Inplace multiply by an arithmetic value
     * 
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return a reference to the updated object
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    MutationalSignatureExprResult& operator*(const T& value)
    {
        if (value>1 || value<0) {
            std::ostringstream oss;

            oss << "the multiplying value must be in the real interval [0,1]. " 
                << std::to_string(value) << " passed.";

            throw std::domain_error(oss.str());
        }

        for (auto& [type, type_value]: value_map) {
            type_value *= value;
        }

        return *this;
    }

    /**
     * @brief Inplace add a mutational signature expression value
     * 
     * @param expression_value is a mutational signature expression value
     * @return a reference to the updated object
     */
    MutationalSignatureExprResult& operator+(MutationalSignatureExprResult&& expression_value);

    /**
     * @brief Inplace add a signature
     * 
     * @param signature is a mutational signature
     * @return a reference to the updated object
     */
    MutationalSignatureExprResult& operator+(const MutationalSignature& signature);

    friend class MutationalSignature;
};

/**
 * @brief A class to represent a mutational signature
 * 
 * A mutational signature is a probability distribution on 
 * the set of mutational types.
 */
class MutationalSignature
{
    std::map<MutationalType, double> dist_map; //!< the signature probability distribution map
public:
    using const_iterator = std::map<MutationalType, double>::const_iterator;

    /**
     * @brief The empty constructor
     */
    MutationalSignature();

    /**
     * @brief A constructor
     * 
     * @param distribution is a mutation type-value map representing a distribution
     */
    explicit MutationalSignature(const std::map<MutationalType, double>& distribution);

    /**
     * @brief Get the initial constant iterator
     * 
     * @return the initial constant iterator 
     */
    inline const_iterator begin() const
    {
        return dist_map.begin();
    }

    /**
     * @brief Get the final constant iterator
     * 
     * @return the final constant iterator 
     */
    inline const_iterator end() const
    {
        return dist_map.end();
    }

    /**
     * @brief Get the probability associated to a mutational type
     * 
     * @param type is the mutational type whose probability is aimed
     * @return the probability of `type`
     */
    double operator()(const MutationalType& type) const;

    /**
     * @brief Multiply by an arithmetic value
     * 
     * @tparam T is the type of the multiplicand
     * @param value is the multiplicand
     * @return the resulting mutational signature expression value
     */
    template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
    inline MutationalSignatureExprResult operator*(const T& value) const
    {
        return MutationalSignatureExprResult(dist_map) * value;
    }

    /**
     * @brief Add a signature
     * 
     * @param signature is a mutational signature
     * @return the resulting mutational signature expression value
     */
    inline MutationalSignatureExprResult operator+(const MutationalSignature& signature) const
    {
        return MutationalSignatureExprResult(dist_map) + signature;
    }

    /**
     * @brief Read mutational signature from a input stream
     * 
     * @param in is the input stream
     * @return a map that associates the name of the signatures in the file 
     *         and the corresponding signature.
     */
    static std::map<std::string, MutationalSignature> read_from_stream(std::istream& in);

    /**
     * @brief Read mutational signature from a input stream
     * 
     * @param in is the input stream
     * @param signature_names is the set of the requested signature
     * @return a map that associates the name of the signatures in the file that
     *         match `signature_names` and the corresponding signature.
     */
    static std::map<std::string, MutationalSignature> read_from_stream(std::istream& in, const std::set<std::string>& signature_names);
};

/**
 * @brief Multiply an arithmetic value and a mutational signature 
 * 
 * @tparam T is the type of the arithmetic value
 * @param value is the arithmetic value
 * @param signature 
 * @return a `MutationalSignatureExprResult` object representing the
 *         the multiplication result 
 */
template<typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
inline MutationalSignatureExprResult operator*(const T& value, const MutationalSignature& signature)
{
    return signature * value;
}

}   // Passengers

}   // Races


namespace std 
{
/**
 * @brief Stream the mutational type in a stream
 * 
 * @param out is the output stream
 * @param type is the mutational type to stream
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Passengers::MutationalType& type);

/**
 * @brief Stream the mutational type from a stream
 * 
 * @param in is the input stream
 * @param type is the object where the streamed mutational type will be placed
 * @return a reference to the input stream
 */
std::istream& operator>>(std::istream& in, Races::Passengers::MutationalType& type);

}  // std

#endif // __RACES_SNV_DISTRIBUTION__
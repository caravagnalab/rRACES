/**
 * @file driver_genotype.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines driver genotype representation
 * @version 0.13
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

#ifndef __RACES_DRIVER_GENOTYPE__
#define __RACES_DRIVER_GENOTYPE__

#include <map>
#include <string>
#include <sstream>

#include "archive.hpp"
#include "cell_event.hpp"

namespace Races 
{

namespace Drivers
{

/**
 * @brief A class to represent promoter methylation/demethylation rates
 */
class EpigeneticRates
{
    double methylation;         //!< the methylation rate
    double demethylation;       //!< the demethylation rate
    
public:
    /**
     * @brief A constructor
     * 
     * @param methylation_rate is the promoter methylation rate
     * @param demethylation_rate is the promoter demethylation rate
     */
    EpigeneticRates(const double methylation_rate, const double demethylation_rate);

    /**
     * @brief A constructor
     * 
     * @param rate is the promoter methylation/demethylation rate
     */
    explicit EpigeneticRates(const double rate);

    /**
     * @brief Get the methylation rate
     * 
     * @return a constant reference to the methylation rate
     */
    inline const double& get_methylation_rate() const
    {
        return methylation;
    }

    /**
     * @brief Set the methylation rate
     * 
     * @param rate is the new methylation rate
     * @return a reference to updated object
     */
    EpigeneticRates& set_methylation_rate(const double& rate);

    /**
     * @brief Get the demethylation rate
     * 
     * @return a constant reference to the demethylation rate
     */
    inline const double& get_demethylation_rate() const
    {
        return demethylation;
    }

    /**
     * @brief Set the demethylation rate
     * 
     * @param rate is the new demethylation rate
     * @return a reference to updated object
     */
    EpigeneticRates& set_demethylation_rate(const double& rate);
};

/**
 * @brief The methylation signature type
 * 
 * The methylation signature represents the epigenetic status 
 * of all the considered promoters. A promoter is set to be 
 * `true` when is methylated and `false` when it is not.
 */
typedef std::vector<bool> MethylationSignature;

class Genotype;

namespace Simulation
{

    class Species;

}   // Simulation

/**
 * @brief A class to represent driver epigenetic genotype
 * 
 * A driver genotype is characterized by specific genomic
 * mutations. The same driver genotype is associated to 
 * different rates depending on its methylation signature.
 * An driver epigenetic genotype is a driver genotype 
 * together with its methylation signature. 
 */
class EpigeneticGenotype {
    static unsigned int counter;                          //!< Total number of driver epigenetic genotype along the computation

    EpigeneticGenotypeId id;                              //!< Identifier
    GenotypeId genomic_id;                                //!< Genomic genotype identifier

    std::string name;                                     //!< Epigenetic genotype name
    MethylationSignature methylation_signature;           //!< Methylation signature

    std::map<CellEventType, double> event_rates;          //!< Event rates

    std::map<EpigeneticGenotypeId, double> epigenetic_rates;  //!< Epigenetic mutation rates 

    /**
     * @brief The empty constructor
     */
    EpigeneticGenotype();

    /**
     * @brief A constructor
     * 
     * @param genotype is the genomic genotype of the new object
     * @param num_of_promoters is the number of methylable promoters
     */
    EpigeneticGenotype(const Genotype& genotype,
                       const size_t num_of_promoters);
public:
    /**
     * @brief Get the driver epigenetic genotype identifier
     * 
     * @return a constant reference to the identifier
     */
    inline const EpigeneticGenotypeId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the driver genomic genotype identifier
     * 
     * @return a constant reference to the genomic genotype identifier
     */
    inline const GenotypeId& get_genomic_id() const
    {
        return genomic_id;
    }

    /**
     * @brief Get the driver epigenetic genotype name
     * 
     * @return a constant reference to the driver genotype name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get the rate of an event 
     * 
     * @param event is the event whose rate is required
     * @return if the rate of `event` has been set, then 
     *       the rate of `event`. The value 0 otherwise
     */
    double get_rate(const CellEventType& event) const;

    /**
     * @brief Set the rate of an event 
     * 
     * @param event is the event whose rate is set
     * @param rate is the new rate for the event
     * @return a constant reference to the new rate
     */
    inline const double& set_rate(const CellEventType& event, const double rate)
    {
        return (event_rates[event] = rate);
    }

    /**
     * @brief Get event rates
     * 
     * @return a constant reference to event rates
     */
    inline const std::map<CellEventType, double>& get_rates() const
    {
        return event_rates;
    }

    /**
     * @brief Set event rates
     * 
     * @param rates are the new event rates
     */
    inline void set_rates(const std::map<CellEventType, double>& rates)
    {
        event_rates = rates;
    }

    /**
     * @brief Get the epigenetic mutation rates
     * 
     * @return a constant reference to a map associating the identifier of 
     *      the epigenetic genotype reachable with an epigenetic mutation and 
     *      its rate
     */
    inline const std::map<EpigeneticGenotypeId, double>& get_epigenetic_mutation_rates() const
    {
        return epigenetic_rates;
    }

    /**
     * @brief Get the methylation signature
     * 
     * @return the methylation signature
     */
    inline const MethylationSignature& get_methylation_signature() const
    {
        return methylation_signature;
    }

    /**
     * @brief Save the epigenetic genotype in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & id 
                & genomic_id 
                & name 
                & methylation_signature 
                & event_rates
                & epigenetic_rates;
    }

    /**
     * @brief Load a epigenetic genotype from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded epigenetic genotype 
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static EpigeneticGenotype load(ARCHIVE& archive)
    {
        EpigeneticGenotype genotype;

        archive & genotype.id 
                & genotype.genomic_id 
                & genotype.name 
                & genotype.methylation_signature 
                & genotype.event_rates
                & genotype.epigenetic_rates;

        if (EpigeneticGenotype::counter < static_cast<unsigned int>(genotype.id+1)) {
            EpigeneticGenotype::counter = genotype.id+1;
        }

        return genotype;
    }

    friend class Genotype;
    friend class Simulation::Species;
};

/**
 * @brief Driver genotype class
 */
class Genotype {
    static unsigned int counter;                        //!< Total number of driver genotype along the computation

    EpigeneticGenotypeId id;                                //!< Identifier
    std::string name;                                   //!< Genotype name

    std::vector<EpigeneticGenotype> e_genotypes;        //!< Epigenetic genotypes

    template<typename SIGNATURE_TYPE>
    void validate_signature(const SIGNATURE_TYPE& signature) const;

public:

    /**
     * @brief A constructor
     * 
     * This constructor builds a driver genotype and its epigenetic genotypes. 
     * The number of epigenetic genotypes is determined by the size of the 
     * epigenetic event rate vector. Each cell of the vector correspond to 
     * a pair methylation/demethylation events for a specific target promoter.
     * Since each target promoter can be either methylated ("+") or 
     * non-methylated ("-"), the number of epigenetic genotypes associated to
     * a driver genotype is exponential in the number of the promoters and, 
     * as a consequence, in the number of elements in the epigenetic event 
     * rate vector.
     * 
     * @param name is the name of the driver genotype
     * @param epigenetic_event_rates is the rates of the methylation/demethylation events 
     *          on each target promoters
     */
    Genotype(const std::string& name, const std::vector<EpigeneticRates>& epigenetic_event_rates);

    /**
     * @brief A constructor
     * 
     * This constructor builds a driver genotype having one epigenetic 
     * genotypes.
     * 
     * @param name is the name of the driver genotype
     */
    explicit Genotype(const std::string& name);

    /**
     * @brief Get the driver epigenetic genotype identifier
     * 
     * @return a constant reference to the identifier
     */
    inline const GenotypeId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the driver genotype name
     * 
     * @return a constant reference to the driver genotype name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get the epigenetic genotypes associated to this genomic genotype
     * 
     * @return a constant reference to the epigenetic genotypes
     */
    inline const std::vector<EpigeneticGenotype>& epigenetic_genotypes() const
    {
        return e_genotypes;
    }

    /**
     * @brief Get the number of methylable promoters in the genotype
     * 
     * @return the number of methylable promoters in the genotype
     */
    size_t num_of_promoters() const;

    /**
     * @brief Get the epigenetic genotype associated to a methylation signature
     * 
     * @param methylation_signature is the methylation signature of the aimed epigenetic genotype
     * @return a non-constant reference to the epigenetic genotype associated 
     *      to `methylation_signature`
     */
    EpigeneticGenotype& operator[](const MethylationSignature& methylation_signature);

    /**
     * @brief Get the epigenetic genotype associated to a methylation signature
     * 
     * @param methylation_signature is the methylation signature of the aimed epigenetic genotype
     * @return a constant reference to the epigenetic genotype associated 
     *      to `methylation_signature`
     */
    const EpigeneticGenotype& operator[](const MethylationSignature& methylation_signature) const;

    /**
     * @brief Get the epigenetic genotype associated to a methylation signature
     * 
     * @param methylation_signature is a string representing the signature of the aimed epigenetic genotype
     * @return a non-constant reference to the epigenetic genotype associated 
     *      to `methylation_signature`
     */
    EpigeneticGenotype& operator[](const std::string& methylation_signature);

    /**
     * @brief Get the epigenetic genotype associated to a methylation signature
     * 
     * @param methylation_signature is a string representing the signature of the aimed epigenetic genotype
     * @return a constant reference to the epigenetic genotype associated 
     *      to `methylation_signature`
     */
    const EpigeneticGenotype& operator[](const std::string& methylation_signature) const;

    static size_t string_to_index(const std::string& methylation_signature);
    static std::string index_to_string(const size_t& index, const size_t num_of_promoters);

    static size_t signature_to_index(const MethylationSignature& methylation_signature);
    static MethylationSignature index_to_signature(const size_t& index, const size_t num_of_promoters);

    inline static std::string signature_to_string(const MethylationSignature& methylation_signature)
    {
        const auto signature_index = signature_to_index(methylation_signature);

        return Genotype::index_to_string(signature_index, methylation_signature.size());
    }
};

template<typename SIGNATURE_TYPE>
void Genotype::validate_signature(const SIGNATURE_TYPE& signature) const
{
    const size_t promoters = num_of_promoters();
    if (signature.size()!=promoters) {
        std::stringstream oss;

        oss << "Wrong signature size: " 
            << get_name() << " has ";

        switch (promoters) {
            case 0:
                oss << "no methylable promoters";
                break;
            case 1:
                oss << "1 methylable promoter";
                break;
            default:
                oss << promoters << " methylable promoters";
        }

        throw std::domain_error(oss.str());
    }
}

}   // Drivers

}   // Races

namespace std
{

/**
 * @brief Write information about an epigenetic rates in an output stream
 * 
 * @param out is the output stream
 * @param epigenetic_rates are the epigenetic rates to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Drivers::EpigeneticRates& epigenetic_rates);

/**
 * @brief Write information about an epigenetic genotype in an output stream
 * 
 * @param out is the output stream
 * @param genotype is the epigenetic genotype to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Drivers::EpigeneticGenotype& genotype);

/**
 * @brief Write information about a genomic genotype in an output stream
 * 
 * @param out is the output stream
 * @param genotype is the genomic genotype to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Drivers::Genotype& genotype);

}   // std


#endif // __RACES_DRIVER_GENOTYPE__
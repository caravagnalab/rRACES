/**
 * @file mutational_properties.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a class to represent the mutational properties
 * @version 0.4
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

#ifndef __RACES_MUTATIONAL_PROPERTIES__
#define __RACES_MUTATIONAL_PROPERTIES__

#include <map>
#include <list>
#include <string>

#include "simulation.hpp"
#include "snv.hpp"
#include "cna.hpp"

namespace Races
{

namespace Passengers
{

/**
 * @brief A class representing the passenger mutational properties of a species
 * 
 */
class SpeciesMutationalProperties
{
    /**
     * @brief The passenger mutational properties of an epigenetic genotype
     */
    struct EpigeneticGenotypeProperties
    {
        std::string name;                       //!< the species name
        double mu;                              //!< the species passenger mutation rate
        std::list<SNV> SNVs;                    //!< the driver species SNVs
        std::list<CopyNumberAlteration> CNAs;   //!< the driver species CNAs

        /**
         * @brief The empty constructor
         */
        EpigeneticGenotypeProperties();

        /**
         * @brief A constructor
         * 
         * @param name is the species name
         * @param SNVs is the vector of species specific SNVs
         * @param CNAs is the vector of species specific CNAs
         */
        EpigeneticGenotypeProperties(const std::string& name, const double& mu,
                                     const std::list<SNV>& SNVs,
                                     const std::list<CopyNumberAlteration>& CNAs);
    };

    std::map<Drivers::EpigeneticGenotypeId, EpigeneticGenotypeProperties> properties;  //!< the epigenetic species property map
public:

    /**
     * @brief The empty constructor
     */
    SpeciesMutationalProperties();

    /**
     * @brief Add new species mutational properties
     * 
     * @param simulation is a driver tissue simulation
     * @param name is the name of the mutational properties genotype
     * @param epigenetic_status_id is a map from epigenomic status to 
     *          passenger mutational rate
     * @param species_SNVs is a list of SNVs characterizing the species
     * @param species_CNAs is a list of CNAs characterizing the species
     * @return a reference to the updated object
     */
    SpeciesMutationalProperties& add_species(const Drivers::Simulation::Simulation& simulation,
                                             const std::string& name, 
                                             const std::map<std::string, double>& epigenetic_status_id,
                                             const std::list<SNV>& species_SNVs={},
                                             const std::list<CopyNumberAlteration>& species_CNAs={});

    /**
     * @brief Get the properties of a species
     * 
     * @param id is the epigenetic identifier of the species
     * @return the mutational properties of the epigenetic species having 
     *          `id` as id.
     */
    inline const SpeciesMutationalProperties::EpigeneticGenotypeProperties& at(const Drivers::EpigeneticGenotypeId& id) const
    {
        return properties.at(id);
    }
};


}   // Passengers

}   // Races

#endif // __RACES_MUTATIONAL_PROPERTIES__
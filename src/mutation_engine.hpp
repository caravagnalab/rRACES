/*
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023-2024 Alberto Casagrande <alberto.casagrande@uniud.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __RRACES_MUTATION_ENGINE__
#define __RRACES_MUTATION_ENGINE__

#include <string>

#include <Rcpp.h>

#include <mutation_engine.hpp>

#include "cna.hpp"
#include "snv.hpp"
#include "phylogenetic_forest.hpp"
#include "samples_forest.hpp"

#include "genomic_data_storage.hpp"

class MutationEngine
{
    using AbsGenotypePosition = uint32_t;

    std::filesystem::path reference_path;

    Races::Mutations::ContextIndex<AbsGenotypePosition> context_index;
    Races::Mutations::MutationEngine<AbsGenotypePosition, std::mt19937_64> m_engine;

    void init_mutation_engine(const GenomicDataStorage& storage,
                              const size_t& default_num_of_alleles,
                              const std::map<std::string, size_t>& alleles_num_exceptions,
                              const size_t& context_sampling_rate=100,
                              const std::string& tumor_type="");
public:
    MutationEngine(const std::string& setup_code,
                   const size_t& context_sampling_rate=100,
                   const std::string& tumor_type="");

    MutationEngine(const std::string& directory,
                   const std::string& reference_url,
                   const std::string& SBS_url,
                   const size_t& default_num_of_alleles,
                   const std::map<std::string, size_t>& exceptions_on_allele_number,
                   const size_t& context_sampling=100,
                   const std::string& tumor_type="");

    static Rcpp::List get_supported_setups();

    static std::list<std::string> get_supported_tumor_types();

    void add_exposure(const Rcpp::List& exposure);

    void add_exposure(const double& time, const Rcpp::List& exposure);

    void add_mutant(const std::string& mutant_name, const Rcpp::List& passenger_rates,
                    const Rcpp::List& driver_SNVs);

    void add_mutant(const std::string& mutant_name, const Rcpp::List& passenger_rates,
                    const Rcpp::List& driver_SNVs, const Rcpp::List& driver_CNAs);

    PhylogeneticForest place_mutations(const SamplesForest& forest,
                                       const size_t& num_of_preneoplatic_mutations,
                                       const int seed);

    inline PhylogeneticForest place_mutations(const SamplesForest& forest,
                                              const size_t& num_of_preneoplatic_mutations)
    {
        return place_mutations(forest, num_of_preneoplatic_mutations, 0);
    }

    inline std::filesystem::path get_reference_path() const
    {
        return reference_path;
    }

    void show() const;

    static MutationEngine 
    build_MutationEngine(const std::string& directory,
                         const std::string& reference_url,
                         const std::string& SBS_url,
                         const size_t& default_num_of_alleles,
                         const Rcpp::List& exceptions_on_allele_number,
                         const std::string& setup_code,
                         const size_t& context_sampling,
                         const std::string& tumor_type);
};

RCPP_EXPOSED_CLASS(MutationEngine)

#endif // __RRACES_SETUP_MUTATION_ENGINE__

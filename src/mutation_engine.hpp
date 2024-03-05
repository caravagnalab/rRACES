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

    GenomicDataStorage storage;

    std::string germline_subject;
    size_t context_sampling;
    std::string tumor_type;

    Races::Mutations::ContextIndex<AbsGenotypePosition> context_index;
    Races::Mutations::MutationEngine<AbsGenotypePosition, std::mt19937_64> m_engine;

    GermlineSubject get_germline_subject(const std::string& subject_name) const;

    void init_mutation_engine();
public:
    MutationEngine(const std::string& setup_code,
                   const std::string& germline_subject="",
                   const size_t& context_sampling=100,
                   const std::string& tumor_type="");

    MutationEngine(const std::string& directory,
                   const std::string& reference_source,
                   const std::string& SBS_source,
                   const std::string& drivers_source,
                   const std::string& passenger_CNAs_source,
                   const std::string& germline_source,
                   const std::string& germline_subject="",
                   const size_t& context_sampling=100,
                   const std::string& tumor_type="");

    static Rcpp::List get_supported_setups();

    void add_exposure(const Rcpp::List& exposure);

    void add_exposure(const double& time, const Rcpp::List& exposure);

    void add_mutant(const std::string& mutant_name, const Rcpp::List& passenger_rates,
                    const Rcpp::List& driver_SNVs);

    void add_mutant(const std::string& mutant_name, const Rcpp::List& passenger_rates,
                    const Rcpp::List& driver_SNVs, const Rcpp::List& driver_CNAs);

    inline Rcpp::List get_active_germline() const
    {
        return storage.get_germline_storage().get_subject_df(germline_subject);
    }

    inline Rcpp::List get_germline_subjects() const
    {
        return storage.get_germline_storage().get_population_df();
    }

    inline Rcpp::List get_population_descritions() const
    {
        return storage.get_germline_storage().get_population_descritions_df();
    }

    inline bool get_infinite_sites_model() const
    {
        return m_engine.infinite_sites_model;
    }

    inline void set_infinite_sites_model(const bool infinite_sites_model)
    {
        m_engine.infinite_sites_model = infinite_sites_model;
    }

    void set_germline_subject(const std::string& germline_subject);

    PhylogeneticForest place_mutations(const SamplesForest& forest,
                                       const size_t& num_of_preneoplatic_mutations,
                                       const int seed);

    inline PhylogeneticForest place_mutations(const SamplesForest& forest,
                                              const size_t& num_of_preneoplatic_mutations)
    {
        return place_mutations(forest, num_of_preneoplatic_mutations, 0);
    }

    Rcpp::List get_SBS_dataframe();

    void show() const;

    static MutationEngine 
    build_MutationEngine(const std::string& directory,
                         const std::string& reference_source,
                         const std::string& SBS_source,
                         const std::string& drivers_source,
                         const std::string& passenger_CNAs_source,
                         const std::string& germline_source,
                         const std::string& setup_code,
                         const std::string& germline_subject,
                         const size_t& context_sampling,
                         const std::string& tumor_type);

    void set_context_sampling(const size_t& context_sampling);

    void rebuild_context_index();

    void reset(const bool full=true);
};

RCPP_EXPOSED_CLASS(MutationEngine)

#endif // __RRACES_SETUP_MUTATION_ENGINE__

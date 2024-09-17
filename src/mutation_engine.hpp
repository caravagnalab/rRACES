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
#include "sid.hpp"
#include "phylogenetic_forest.hpp"
#include "samples_forest.hpp"

#include "genomic_data_storage.hpp"
#include "utility.hpp"

class MutationEngine
{
    using AbsGenotypePosition = uint32_t;

    GenomicDataStorage storage;

    std::string germline_subject;
    size_t context_sampling;
    size_t max_motif_size;
    size_t max_repetition_storage;
    std::string tumour_type;
    std::string tumor_study;

    RACES::Mutations::ContextIndex<AbsGenotypePosition> context_index;
    RACES::Mutations::RSIndex rs_index;
    RACES::Mutations::MutationEngine<AbsGenotypePosition, std::mt19937_64> m_engine;

    bool avoid_homozygous_losses;

    GermlineSubject get_germline_subject(const std::string& subject_name) const;

    void init_mutation_engine(const bool& quiet);

    template<typename MUTATION_TYPE>
    void show_timed_exposures() const
    {
        Rcpp::Rcout << "   " << MUTATION_TYPE::name() << " Timed Exposures" << std::endl;

        const auto& timed_exposures = m_engine.get_timed_exposures(MUTATION_TYPE::type());
        auto coeffs_it = timed_exposures.begin();
        while (coeffs_it != timed_exposures.end()) {
            auto next_it = coeffs_it;
            ++next_it;
            if (next_it == timed_exposures.end()) {
                Rcpp::Rcout << "     [" << coeffs_it->first << ", \u221E[: ";
            } else {
                Rcpp::Rcout << "     ["
                            << coeffs_it->first << ", " << next_it->first << "[: ";
            }
            Rcpp::Rcout << coeffs_it->second << std::endl;

            coeffs_it = next_it;
        }
    }
public:
    MutationEngine(const std::string& setup_code,
                   const std::string& germline_subject="",
                   const size_t& context_sampling=100,
                   const size_t& max_motif_size=50,
                   const size_t& max_repetition_storage=500000,
                   const std::string& tumour_type="",
                   const std::string& tumor_study="",
                   const bool& avoid_homozygous_losses=true,
                   const bool& quiet=false);

    MutationEngine(const std::string& directory,
                   const std::string& reference_source,
                   const std::string& SBS_signatures_source,
                   const std::string& indel_signatures_source,
                   const std::string& drivers_source,
                   const std::string& passenger_CNAs_source,
                   const std::string& germline_source,
                   const std::string& germline_subject="",
                   const size_t& context_sampling=100,
                   const size_t& max_motif_size=50,
                   const size_t& max_repetition_storage=500000,
                   const std::string& tumour_type="",
                   const std::string& tumor_study="",
                   const bool& avoid_homozygous_losses=true,
                   const bool& quiet=false);

    static Rcpp::List get_supported_setups();

    inline void add_exposure(const Rcpp::List& exposure)
    {
        add_exposure(0, exposure);
    }

    void add_exposure(const double& time, const Rcpp::List& exposure);

    void add_mutant(const std::string& mutant_name, const Rcpp::List& passenger_rates);

    void add_mutant(const std::string& mutant_name, const Rcpp::List& passenger_rates,
                    const Rcpp::List& drivers);

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

    Rcpp::List get_known_driver_mutations() const;

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
                                       const size_t& num_of_preneoplatic_SNVs,
                                       const std::string& preneoplatic_SNV_signature_name,
                                       const size_t& num_of_preneoplatic_indels,
                                       const std::string& preneoplatic_indel_signature_name,
                                       const int seed);

    inline PhylogeneticForest place_mutations(const SamplesForest& forest,
                                              const size_t& num_of_preneoplatic_SNVs,
                                              const std::string& preneoplatic_SNV_signature_name,
                                              const size_t& num_of_preneoplatic_indels,
                                              const std::string& preneoplatic_indel_signature_name,
                                              const SEXP& seed)
    {
        return place_mutations(forest, num_of_preneoplatic_SNVs, preneoplatic_SNV_signature_name,
                               num_of_preneoplatic_indels, preneoplatic_indel_signature_name,
                               get_random_seed<int>(seed));
    }

    inline PhylogeneticForest place_mutations(const SamplesForest& forest,
                                              const size_t& num_of_preneoplatic_SNVs,
                                              const std::string& preneoplatic_SNV_signature_name,
                                              const size_t& num_of_preneoplatic_indels,
                                              const std::string& preneoplatic_indel_signature_name)
    {
        return place_mutations(forest, num_of_preneoplatic_SNVs, preneoplatic_SNV_signature_name,
                               num_of_preneoplatic_indels, preneoplatic_indel_signature_name,
                               get_random_seed<int>(R_NilValue));
    }

    inline PhylogeneticForest place_mutations(const SamplesForest& forest,
                                              const size_t& num_of_preneoplatic_SNVs,
                                              const size_t& num_of_preneoplatic_indels,
                                              const SEXP& seed)
    {
        return place_mutations(forest, num_of_preneoplatic_SNVs, "SBS1",
                               num_of_preneoplatic_indels, "ID1", seed);
    }

    inline PhylogeneticForest place_mutations(const SamplesForest& forest,
                                              const size_t& num_of_preneoplatic_SNVs,
                                              const size_t& num_of_preneoplatic_indels)
    {
        return place_mutations(forest, num_of_preneoplatic_SNVs, "SBS1",
                               num_of_preneoplatic_indels, "ID1",
                               get_random_seed<int>(R_NilValue));
    }

    Rcpp::List get_SNV_signatures_dataframe() const;

    Rcpp::List get_indel_signatures_dataframe() const;

    void show() const;

    static MutationEngine
    build_MutationEngine(const std::string& directory,
                         const std::string& reference_source,
                         const std::string& SBS_signatures_source,
                         const std::string& indel_signatures_source,
                         const std::string& drivers_source,
                         const std::string& passenger_CNAs_source,
                         const std::string& germline_source,
                         const std::string& setup_code,
                         const std::string& germline_subject,
                         const size_t& context_sampling,
                         const size_t& max_motif_size,
                         const size_t& max_repetition_storage,
                         const std::string& tumour_type,
                         const std::string& tumour_study,
                         const bool avoid_homozygous_losses,
                         const bool quiet);

    static Rcpp::List get_available_tumour_type(const std::string& setup_code);

    inline void set_context_sampling(const size_t& context_sampling)
    {
        set_context_sampling(context_sampling, false); 
    }

    void set_context_sampling(const size_t& context_sampling,
                              const bool quiet);

    inline void rebuild_indices()
    {
        rebuild_indices(false);
    }

    void rebuild_indices(const bool quiet);

    void reset(const bool full=true, const bool quiet=false);
};

RCPP_EXPOSED_CLASS(MutationEngine)

#endif // __RRACES_SETUP_MUTATION_ENGINE__

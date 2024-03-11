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

#include <string>

#include <Rcpp.h>

#include "phylogenetic_forest.hpp"
#include "mutation_engine.hpp"
#include "snv.hpp"
#include "cna.hpp"

using namespace Rcpp;

RCPP_MODULE(Mutations){

//' @name SNV
//' @title A single nucleotide variation
  class_<SNV>("SNV")
    .constructor()

//' @name SNV$get_chromosome
//' @title Get the chromosome in which the SNV occurs.
//' @return The chromosome in which the SNV occurs.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the chromosome in which `snv` occurs (i.e., "X")
//' snv$get_chromosome()
    .method("get_chromosome",&SNV::get_chromosome, "Get the chromosome of the SNV position")

//' @name SNV$get_position_in_chromosome
//' @title Get the position in chromosome where the SNV occurs.
//' @return The position in chromosome where the SNV occurs.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the position in chromosome where `snv` occurs (i.e., 20002)
//' snv$get_position_in_chromosome()
    .method("get_position_in_chromosome",&SNV::get_position_in_chromosome,
            "Get the SNV position in the chromosome")

//' @name SNV$get_ref_base
//' @title Get the reference base in which the SNV occurs.
//' @return The reference base in which the SNV occurs.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the reference base in which `snv` occurs (i.e., "A")
//' snv$get_ref_base()
    .method("get_ref_base",&SNV::get_ref_base, "Get the reference base in the SNV locus")

//' @name SNV$get_alt_base
//' @title Get the base after the SNV occurs.
//' @return The base after the SNV occurs.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the base after `snv` occurs (i.e., "T")
//' snv$get_alt_base()
    .method("get_alt_base",&SNV::get_alt_base, "Get the new base of the SNV")

//' @name SNV$get_cause
//' @title Get the SNV cause.
//' @description Evey SNV produced by RACES/rRACES is associated to a cause depending
//'   on whether it is part of a genomic characterization of a mutant or it is 
//'   caused by a specific SBS profile. This method return such a cause whenever 
//'   it is available.
//' @return The base after the SNV occurs.
//' @examples
//' # let us build a SNV without specifying any cause for it
//' snv <- SNV("X", 20002, "T", "A")
//'
//' # get the cause of `snv` (i.e., "NA")
//' snv$get_cause()
//'
//' # we can also build a SNV, specifying a cause for it
//' snv <- SNV("X", 20002, "T", "A", cause = "SBS13")
//'
//' # get the cause of `snv` (i.e., "SBS13")
//' snv$get_cause()
    .method("get_cause",&SNV::get_cause, "Get the cause of the SNV")

//' @name SNV$get_dataframe
//' @title Get a data frame representing the SNV.
//' @description This method returns a data frame that represents the SNV and whose
//'   columns are `chromosome`, `pos_in_chr`, `ref`, `alt`, and
//'   `cause`.
//' @examples
//' snv <- SNV("X", 20002, "T", "A")
//'
//' snv$get_dataframe()
    .method("get_dataframe",&SNV::get_dataframe, "Get a dataframe representing the SNV")
    .method("show",&SNV::show);

//' @name SNV
//' @title Create a SNV amplification.
//' @param chromosome The name of the chromosome in which the SNV occurs.
//' @param pos_in_chr The position in the chromosome where the SNV occurs.
//' @param alt The base after the mutation.
//' @param ref The base before the mutation (optional).
//' @param cause The cause of the SNV (optional).
//' @examples
//' # create a SNV without specifying the cause and context 
//' snv <- SNV("X", 20002, "T")
//' snv
//'
//' # create a SNV and do not specify the cause
//' snv <- SNV("X", 20002, "T", "A")
//' snv
//'
//' # create a SNV with a cause
//' snv <- SNV("X", 20002, "T", cause = "SBS1")
//' snv
  function("SNV", &SNV::build_SNV,
           List::create(_["chromosome"], _["pos_in_chr"], _["alt"], 
                        _["ref"] = "?", _["cause"] = ""),
           "Create a single nucleotide variation (SNV)");

//' @name CNA
//' @title Create a CNA amplification.
//' @param type The CNA type: either `"A"` or `"D"` for amplification and 
//'   deletion, respectively.
//' @param chromosome The name of the chromosome in which the CNA occurs.
//' @param pos_in_chr The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the CNA occurs. (optional)
//' @param src_allele The allele from which the region is amplified. (optional, 
//'   for amplification only)
//' @seealso `Amplification` to build an amplification; `Deletion` to build a 
//'   deletion.
//' @examples
//' # create an amplification
//' cna <- CNA("A", "X", 20002, 100)
//'
//' cna
//'
//' # create a deletion from the allele 0
//' cna <- CNA("D", "Y", 101310, 205, allele=0)
//' 
//' cna
  function("CNA", &CNA::build_CNA,
           List::create(_["type"], _["chromosome"], _["pos_in_chr"], _["len"],
                        _["allele"]=false, _["src_allele"]=false),
           "Create a copy number alteration (CNA)");

//' @name Amplification
//' @title Create a CNA amplification.
//' @param chromosome The name of the chromosome in which the CNA occurs.
//' @param pos_in_chr The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the amplification is placed. (optional)
//' @param src_allele The allele from which the region is amplified. (optional)
//' @seealso `Deletion` to build a deletion; `CNA` to build both amplifications 
//'   and deletions.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100)
//'
//' cna
  function("Amplification", &CNA::build_amplification,
           List::create(_["chromosome"], _["pos_in_chr"], _["len"],
                        _["allele"]=false, _["src_allele"]=false),
           "Create a CNA amplification");

//' @name Deletion
//' @title Create a CNA deletion.
//' @param chromosome The name of the chromosome in which the CNA occurs.
//' @param pos_in_chr The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the deletion occurs. (optional)
//' @seealso `Amplification` to build an amplification; `CNA` to build 
//'   both amplifications and deletions.
//' @examples
//' # create a deletion CNA
//' cna <- Deletion("Y", 40020, 200)
//'
//' cna
  function("Deletion", &CNA::build_deletion,
           List::create(_["chromosome"], _["pos_in_chr"], _["len"],
                        _["allele"]=false),
           "Create a CNA deletion");

//' @name CNA
//' @title A copy number alteration
  class_<CNA>("CNA")
    .constructor()

//' @name CNA$get_chromosome
//' @title Get the chromosome in which the CNA occurs.
//' @return The chromosome in which the CNA occurs.
//' @examples
//' # create an amplification CNA
//' cna <- CNA("A", "X", 20002, 100)
//'
//' # get the chromosome in which `cna` occurs (i.e., "X")
//' cna$get_chromosome()
    .method("get_chromosome",&CNA::get_chromosome, "Get the chromosome of the CNA")

//' @name CNA$get_position_in_chromosome
//' @title Get the position in chromosome where the CNA occurs.
//' @return The position in chromosome where the CNA occurs.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get the position in chromosome where `cna` occurs (i.e., 20002)
//' cna$get_position_in_chromosome()
    .method("get_position_in_chromosome",&CNA::get_position_in_chromosome,
            "Get the CNA position in the chromosome")

//' @name CNA$get_length
//' @title Get the CNA length.
//' @return The CNA length.
//' @examples
//' # create an amplification CNA
//' cna <- CNA("A", "X", 20002, 100)
//'
//' # get the length of `cna` (i.e., 100)
//' cna$get_length()
    .method("get_length",&CNA::get_length, "Get the CNA length")

//' @name CNA$get_allele
//' @title Get the allele in which CNA occurs.
//' @return The allele in which CNA occurs.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get the allele in which `cna` occurs (i.e., 1)
//' cna$get_allele()
    .method("get_allele",&CNA::get_allele, "Get the alteration allele")

//' @name CNA$get_src_allele
//' @title Get the allele from which CNA is copied.
//' @return The allele from which CNA is copied.
//' @examples
//' # create an amplification CNA
//' amp_cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' # get allele from which `amp_cna` is copied (i.e., 0)
//' amp_cna$get_src_allele()
//'
//' # create a deletion CNA
//' del_cna <- Deletion("Y", 40020, 200, 0)
//'
//' # the deletions have no sources and the method returns NA
//' del_cna$get_src_allele()
    .method("get_src_allele",&CNA::get_src_allele, "Get the source allele (for amplifications)")

//' @name CNA$get_dataframe
//' @title Get a data frame representing the CNA.
//' @description This method returns a data frame that represents the CNA and whose
//'   columns are `chromosome`, `pos_in_chr`, `length`, `alt_base`, `allele`,
//'   `src_allele`, and `type`.
//' @examples
//' # create an amplification CNA
//' amp_cna <- Amplification("X", 20002, 100)
//'
//' amp_cna$get_dataframe()
//'
//' # create a deletion CNA
//' del_cna <- Deletion("Y", 40020, 200, 0)
//'
//' del_cna$get_dataframe()
    .method("get_dataframe",&CNA::get_dataframe, "Get a dataframe representing the CNA")
    .method("show",&CNA::show);

//' @name MutationEngine
//' @title Mutation engines generate phylogenetic forests
//' @description A mutation engine can label every node of a descendants 
//'   forest by mutations and produce a consistent phylogenetic forest.
//'   The mutations are randomly generated according to three factors:
//'     - the mutational rates of the species involved in the descendants 
//'       forest 
//'     - the genotycal characterization of the mutants involved in the 
//'       descendants forest, i.e., the SVNs and CNAs characterizing the
//'       mutant genotypes
//'     - the SBS signature coefficients active along the species
//'       simulation
//'
//'   These data are provided to a mutation engine by using the methods 
//'   [MutationEngine$add_exposure()] and [MutationEngine$add_exposure()]
//'   These data are provided by means of the [MutationEngine$add_mutant()].
//'   
//'   The construction of a `MutationEngine` object requires a reference 
//'   sequence and an SBS file which are downloaded from the Internet. 
//'   After the download a context index of the reference sequence is then 
//'   automatically built. Thess processes may take time depending on the 
//'   size of the reference sequence. Because of this, the downloaded files 
//'   together with the context index are saved in a directory on the disk 
//'   and they are available for successive `MutationEngine` constructions.
//'   
  class_<MutationEngine>("MutationEngine")

//' @name MutationEngine$infinite_sites_model
//' @title A flag to enable/disable the infinite sites model.
//' @description This property enables/disables the infinite sites model.
//'   When it is `TRUE`, the infinite sites model is enabled and 
//'   new mutations are exclusively placed in mutation-free loci.
//' @examples
//' # create a demostrative mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # the infinite sites model is enabled by default
//' m_engine$infinite_sites_model
//' 
//' # the infinite sites model can be disabled
//' m_engine$infinite_sites_model <- FALSE
//'
//' m_engine$infinite_sites_model
    .property("infinite_sites_model", &MutationEngine::get_infinite_sites_model, 
              &MutationEngine::set_infinite_sites_model,
              "A flag to enable/disable the infinite sites model")

//' @name MutationEngine$add_exposure
//' @title Add an exposure to the mutation engine. 
//' @description This method adds an exposure to the mutation engine.
//'   
//'   The exposure will be used to establish the probability 
//'   for a passenger mutation to occur depending on its context.
//'
//'   Each exposure is associated to a time that is the simulated
//'   time in which the set is adopted.
//'   If a time is provided the exposure is used from the
//'   specified time on up to the successive possible change of
//'   exposure. When an exposure is added to the mutation engine
//'   without specifying the time, its time is 0.
//' @param time The simulated time at which the exposure is adopted 
//'   (optional).
//' @param exposure An signature exposure, i.e., a discrete probability
//'   distribution over a set of SBS signature.
//' @examples
//' # create a demostrative mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # add a default set of SBS coefficients that will be used from simulated
//' # time 0 up to the successive SBS coefficient change.
//' m_engine$add_exposure(c(SBS13 = 0.3, SBS1 = 0.7))
//'
//' # add a default set of SBS coefficients that will be used from simulated
//' # time 3.2 up to the end of the simulation.
//' m_engine$add_exposure(3.2, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5))
//'
//' m_engine
    .method("add_exposure", (void (MutationEngine::*)(const List&))(
              &MutationEngine::add_exposure), "Add an exposure")
    .method("add_exposure", (void (MutationEngine::*)(const double&, const List&))(
              &MutationEngine::add_exposure), "Add an exposure")

//' @name MutationEngine$add_mutant
//' @title Add a mutant specification to the mutation engine.
//' @description This method adds a mutant specification to the mutation engine.
//'   The users must use it to specify the name and the genomic
//'   characterization (i.e., SNVs and CNAs) of all the simulated mutants
//'   together with the mutation rates of its species.
//' @param mutant_name The mutant name.
//' @param passenger_rates The list of the passenger rates whose names are the 
//'   epigenetic states of the species or a single rate, if the mutant
//'   does not have an epigenetic state.
//' @param driver_SNVs The list of the driver SNVs characterizing the mutant.
//' @param driver_CNAs The list of the driver CNAs characterizing the mutant.
//' @examples
//' # create a demostrative mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # add the mutant "A" characterized by one driver SNV on chromosome 22 and
//' # two CNAs: an amplification and a deletion. The mutant has two epigenetic
//' # states and its species "A+" and "A-" have passenger SNV rates 1e-9 and
//' # 3e-8, respectively, and passenger CNA rates 0 and 1e-11, respectively.
//' m_engine$add_mutant("A", list("+" = c(SNV = 1e-9),
//'                               "-" = c(SNV = 3e-8, CNA = 1e-11)),
//'                     driver_SNVs = c(SNV("22", 10510210, "C")),
//'                     driver_CNAs = c(CNA(type = "A", "22", 
//'                                         pos_in_chr = 10303470,
//'                                         len = 200000),
//'                                     CNA("D", "22", 5010000, 200000)))
//'
//' # add the mutant "B" characterized by one driver SNV on chromosome 1 (no
//' # CNA) and missing epigenetic state. Its species "B" has passenger SNV
//' # rate 5e-9 and passenger CNA rate 0.
//' m_engine$add_mutant("B", c(SNV = 5e-9), c(SNV("22", 10510210, "C")))
//'
//' m_engine
    .method("add_mutant", (void (MutationEngine::*)(const std::string&, const Rcpp::List& passenger_rates,
                                                    const Rcpp::List&))(
                                                        &MutationEngine::add_mutant),                                            
            "Add mutant")
    .method("add_mutant", (void (MutationEngine::*)(const std::string&, const Rcpp::List& passenger_rates,
                                                    const Rcpp::List&, const Rcpp::List&))(
                                                        &MutationEngine::add_mutant), 
            "Add mutant")

//' @name MutationEngine$place_mutations
//' @title Place the mutations on a samples forest
//' @description This method labels each node of a samples forest
//'   by the mutations occuring for the first time in the
//'   cell represented by the node itself and produces a
//'   phylogenetic forest.
//' @param samples_forest A samples forest.
//' @param num_of_preneoplatic_mutations The number of pre-neoplastic
//'   mutations.
//' @param seed The seed for random number generator. (optional)
//' @return A phylogenetic forest whose structure corresponds to
//'   `samples_forest`.
//' @examples
//' # create a simulation
//' sim <- new(Simulation)
//' sim$add_mutant("A", c(SNV = 0.2), 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner = c(450, 475),
//'                        upper_corner = c(500, 550))
//'
//' # build the samples forest
//' samples_forest <- sim$get_samples_forest()
//'
//' # build a mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # add the mutant "A" to the engine
//' m_engine$add_mutant("A", c(SNV = 3e-9),
//'                     c(SNV("22", 12028576, "G")))
//'
//' # add the default set of SBS coefficients
//' m_engine$add_exposure(c(SBS13 = 0.3, SBS1 = 0.7))
//'
//' # place the mutations on the samples forest assuming 
//' # 1000 pre-neoplastic mutations
//' phylogenetic_forest <- m_engine$place_mutations(samples_forest, 1000)
//'
//' phylogenetic_forest
    .method("place_mutations", (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest,
                                                                       const size_t& num_of_preneoplatic_mutations))(
                                                        &MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")
    .method("place_mutations", (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest,
                                                                       const size_t& num_of_preneoplatic_mutations,
                                                                       const int seed))(
                                                        &MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")

//' @name MutationEngine$get_active_germline
//' @title Get a data frame describing the active germline subject
//' @description This method returns a data frame containing the
//'   active germline subject. The column `sample` reports the
//'   subject name; the columns `pop` and `super_pop` contain the 
//'   subject population and super population, respectively;
//'   the column `gender` declares the subject gender.
//' @return A data frame the active germline subject.
//' @seealso [MutationEngine$get_germline_subjects()] to get the
//'   available germline subjects; [MutationEngine$set_germline_subject()]
//'   to set the active germline subject.
//' @examples
//' # build a mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_active_germline(), 5)
    .method("get_active_germline", &MutationEngine::get_active_germline)

//' @name MutationEngine$set_germline_subject
//' @title Set the germline subject
//' @description This method sets the germline subject. The subject 
//'   must be one among those reported by 
//'   [MutationEngine$get_germline_subjects()].
//' @return Set the germline subject.
//' @seealso [MutationEngine$get_germline_subjects()] to get the
//'   available germline subjects; [MutationEngine$get_active_germline()]
//'   to get the active germline subject.
//' @examples
//' # build a mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # set the active germline subject data frame
//' m_engine$set_germline_subject("NA18941")
    .method("set_germline_subject", &MutationEngine::set_germline_subject)

//' @name MutationEngine$get_germline_subjects
//' @title Get a data frame reporting the available germline subjects
//' @description This method returns a data frame containing the
//'   avaiable germline subjects. The column `sample` reports the
//'   subject name; the columns `pop` and `super_pop` contain the 
//'   subject population and super population, respectively;
//'   the column `gender` declares the subject gender.
//' @return A data frame the avaiable germline subjects.
//' @seealso [MutationEngine$get_active_germline()] to get the
//'   available germline subjects; [MutationEngine$set_germline_subject()]
//'   to set the active germline.
//' @examples
//' # build a mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_germline_subjects(), 5)
    .method("get_germline_subjects", &MutationEngine::get_germline_subjects)

//' @name MutationEngine$get_population_descritions
//' @title Get a data frame containing the population descriptions
//' @description This method returns a data frame describing the
//'   populations. The column `code` contains the population codes;
//'   the columns `description` and `long description` report a
//'   brief and a long description for the populations,
//'   respectively.
//' @return A data frame containing the population descriptions.
//' @examples
//' # build a mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # get the active germline subject data frame
//' head(m_engine$get_population_descritions(), 5)
    .method("get_population_descritions", &MutationEngine::get_population_descritions)

//' @name MutationEngine$get_SBSs
//' @title Get a data frame containing the available SBSs
//' @description This method returns a data frame containing the
//'   available SBSs and the corresponding mutation
//'   probability. The first column ("Type") describes a 
//'   mutation in a context, while each of the remaining
//'   columns contains the probabilities of the mutations
//'   for one of the available SBSs.
//' @return A data frame containing the available SBSs.
//' @examples
//' # build a mutation engine
//' m_engine <- build_mutation_engine(setup_code = "demo")
//'
//' # get the SBS data frame
//' head(m_engine$get_SBSs(), 5)
    .method("get_SBSs", &MutationEngine::get_SBS_dataframe,
            "Get the SBS data frame")

    .method("show", &MutationEngine::show);

//' @name build_mutation_engine
//' @title Create a mutation engine object
//' @details This function downloads and sets up the data requires by
//'   by a mutation engine object and build it.
//'   
//'   There are two building modalities: the first one is more general, 
//'   but it requires to specify all the data sources; the second one  
//'   adopts some pre-set configurations, but it is more convient than
//'   the former in many cases.
//'
//'   The first building modality requires to specify the directory in
//'   which the data must be saved, the path or URL of the reference 
//'   sequence, the SBS file, the driver SNVs file, the passenger CNAs
//'   file, and the germline data directory through the parameters 
//'   `directory`, `reference_src`, `SBS_src`, `drivers_src`, 
//'   `passenger_CNAs_src`, and `germline_src`, respectively.
//'
//'   The second building modality exclusively requires a set-up code 
//'   (parameter `setup_code`). The list of supported set-up codes can 
//'   be obtained by using the function [get_mutation_engine_codes()].
//'
//'   The number of context sampling is an optional parameter that allows
//'   sampling the reference contexts while building the context index.
//'   This parameter, which is set to 100 by default, specifies how many
//'   occurences of the same context must be identified before adding 
//'   one of them to the context index. The larger the number of context
//'   sampling, the larger the context index. On the other side, the
//'   lower the number of context sampling, the lower the number of sites
//'   in the reference genome that can be affected by simulated
//'   mutations.
//'
//'   If the parameters of a mutation engine construction match those 
//'   of a previous construction, then the corresponding reference
//'   sequence, the SBS file, and the previously built context index
//'   are loaded from the set-up directory avoiding further 
//'   computations.
//' @seealso [get_mutation_engine_codes()] provides a list of the supported 
//'   set-up codes.
//' @export
//' @param setup_code The set-up code (alternative to `directory`).
//' @param directory The set-up directory (alternative to `setup_code`).
//' @param reference_src The reference genome path or URL (mandatory when `directory`
//'   is provided).
//' @param SBS_src The SBS file path or URL (mandatory when `directory` is provided).
//' @param drivers_src The driver mutation file path or URL (mandatory when 
//'   `directory` is provided).
//' @param passenger_CNAs_src The passenger CNAs file path or URL (mandatory when 
//'   `directory` is provided).
//' @param germline_src The germline directory path or URL (mandatory when 
//'   `directory` is provided).
//' @param germline_subject The germline subject (optional).
//' @param context_sampling The number of reference contexts per context in the
//'   index (optional: default value is 100).
//' @param tumor_type The type of tumor. This is currently used to select the
//'   admissible passenger CNAs. If any passenger CNA in the dataset is
//'   admissible, use the the empty string `""` (optional: default value is
//'   `""`).
//' @seealso [MutationEngine$get_germline_subjects()] to get the available
//'   germline subjects; [MutationEngine$set_germline_subject()] to set the
//'   active germline subject; [MutationEngine$get_active_germline()] to get
//'   the active germline subject.
//' @return A mutation engine object.
//' @examples
//' # set the reference and SBS URLs
//' reference_url <- paste0("https://ftp.ensembl.org/pub/grch37/current/",
//'                         "fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.",
//'                         "dna.chromosome.22.fa.gz")
//' sbs_url <- paste0("https://cancer.sanger.ac.uk/signatures/documents/2123/",
//'                   "COSMIC_v3.4_SBS_GRCh37.txt")
//' drivers_url <- paste0("https://raw.githubusercontent.com/",
//'                       "caravagnalab/rRACES/main/inst/extdata/",
//'                       "driver_mutations_hg19.csv")
//' passenger_CNAs_url <- paste0("https://raw.githubusercontent.com/",
//'                              "caravagnalab/rRACES/main/inst/extdata/",
//'                              "passenger_CNAs_hg19.csv")
//' germline_url <- paste0("https://www.dropbox.com/scl/fi/g9oloxkip18tr1r",
//'                        "m6wjve/germline_data_demo.tar.gz?rlkey=15jshul",
//'                        "d3bqgyfcs7fa0bzqeo&dl=1")
//'
//' # build a mutation engine and save the required files into the
//' # directory "Test". The `drivers_url` parameter is optional, but 
//' # it is suggested to avoid passenger mutations on driver loci.
//' m_engine <- build_mutation_engine(directory = "Test",
//'                                   reference_src = reference_url,
//'                                   SBS_src = sbs_url,
//'                                   drivers_src = drivers_url,
//'                                   passenger_CNAs_src = passenger_CNAs_url,
//'                                   germline_src = germline_url)
//'
//' # if the parameters of a mutation engine construction match those of a
//' # previous construction, then the corresponding reference sequence,
//' # the SBS file, and the previously built context index are loaded from
//' # the set-up directory avoiding further computations.
//' m_engine <- build_mutation_engine("Test", reference_url, sbs_url,
//'                                   drivers_url, passenger_CNAs_url,
//'                                   germline_url)
//'
//' # if the `context_sampling` parameter changes, a new context index is
//' # built, while neither the reference sequence nor the SBS file are 
//' # downloaded again.
//' m_engine <- build_mutation_engine("Test", reference_url, sbs_url,
//'                                   drivers_url, passenger_CNAs_url,
//'                                   germline_url, context_sampling = 50)
//'
//' # a futher contruction with the same parameters avoids both downloads
//' # and context index construction.
//' m_engine <- build_mutation_engine("Test", reference_url, sbs_url,
//'                                   drivers_url, passenger_CNAs_url,
//'                                   germline_url, context_sampling = 50)
//'
//' m_engine
//'
//' # the parameters `directory`, `reference_src`, `SBS_src`, `drivers_src`,
//' # `passenger_CNAs_src`, and `germline_src` can be avoided by providing
//' # the `setup_code` parameter. The set-up code `demo` is provided among 
//' # those available for testing purpose.
//' m_engine <- build_mutation_engine(setup_code = "demo")
//' 
//' # the `context_sampling` can be used also when a pre-defined set-up
//' # configuration is adopted.
//' m_engine <- build_mutation_engine(setup_code = "demo",
//'                                   context_sampling = 50)
//'
//' m_engine
//' 
//' # remove the "Test" and "demo" directories
//' unlink("Test", recursive = TRUE)
//' unlink("demo", recursive = TRUE)
  function("build_mutation_engine", &MutationEngine::build_MutationEngine,
           List::create(_["directory"] = "",
                        _["reference_src"] = "", _["SBS_src"] = "",
                        _["drivers_src"] = "", _["passenger_CNAs_src"] = "",
                        _["germline_src"] = "", _["setup_code"] = "", 
                        _["germline_subject"] = "", _["context_sampling"] = 100, 
                        _["tumor_type"] = ""),
           "Create a MutationEngine");

//' @name get_mutation_engine_codes
//' @title Get the supported codes for predefined set-up
//' @return A data frame reporting the code and a description for each
//'   supported predefined set-up.
//' @seealso [build_mutation_engine()] to build a mutation engine
//' @export
//' @examples
//' # get the list of supported mutation engine set-up codes
//' get_mutation_engine_codes()
  function("get_mutation_engine_codes", &MutationEngine::get_supported_setups,
           "Get mutation engine supported codes");

//' @name PhylogeneticForest
//' @title The phylogenetic forest of cells in samples.
//' @description Represents the phylogenetic forest of the
//'   cells sampled during the computation. The leaves of
//'   this forest are the sampled cells.
//'   This class is analoguous to the class `SamplesForest`, 
//'   but each node is labelled with the mutations occuring
//'   for the first time on the cell represented by the node
//'   itself. Moreover each leaf is also associated with the
//'   genome mutations occurring in the corresponding cell.
//' @field get_coalescent_cells Retrieve most recent common ancestors\itemize{
//' \item \emph{Parameter:} \code{cell_ids} - The list of the identifiers of the
//'   cells whose most recent common ancestors are aimed (optional).
//' \item \emph{Return:} A data frame representing, for each of the identified
//'   cells, the identified (column `cell_id`), whenever the
//'   node is not a root, the ancestor identifier (column
//'   `ancestor`), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column `sample`), the mutant
//'   (column `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' }
//' @field get_first_occurrences Gets the identifier of the cell in which a 
//'   mutation occurs for the first time\itemize{
//' \item \emph{Parameter:} \code{mutation} - A mutation being either a 
//'   SNV or a CNA.
//' \item \emph{Return:} The identifier of the cell in which a mutation 
//'   occurs for the first time.
//' }
//' @field get_germline_SNVs Gets the SNVs of the germline\itemize{
//' \item \emph{Return:} A data frame reporting `chromosome`, `chr_pos` (i.e., 
//'   the position in the chromosome), `allele` (in which the SNV
//'   occurs), `ref`, `alt`, and `cause`.
//' }
//' @field get_germline_subject Gets the germline subject name\itemize{
//' \item \emph{Return:} The name of the subject whose germline is used.
//' }
//' @field get_nodes Get the forest nodes \itemize{
//' \item \emph{Return:} A data frame representing, for each node
//'   in the forest, the identified (column `id`),
//'   whenever the node is not a root, the ancestor
//'   identifier (column `ancestor`), whenever the node
//'   was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column `sample`), the mutant (column
//'   `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' }
//' @field get_sampled_cell_CNAs Gets the CNAs of the sampled cells \itemize{
//' \item \emph{Returns:} A data frame reporting `cell_id`, `type` (`"A"` for 
//'   amplifications and `"D"` for deletions), `chromosome`, `begin`
//'   (i.e., the first CNA locus in the chromosome), `end` (i.e., the
//'   last CNA locus in the chromosome), `allele`, and `src allele` 
//'   (the allele origin for amplifications, `NA` for deletions).
//' }
//' @field get_sampled_cell_SNVs Gets the SNVs of the sampled cells \itemize{
//' \item \emph{Returns:} A data frame reporting `cell_id`, `chromosome`, 
//'   `chr_pos` (i.e., position in the chromosome), `allele` (in which
//'   the SNV occurs), `ref`, `alt`, and `cause` for each
//'   SNV in the sampled cell genomes.
//' }
//' @field get_samples_info Retrieve information about the samples \itemize{
//' \item \emph{Returns:} A data frame containing, for each sample collected
//'   during the simulation, the columns `name`, `time`, `ymin`,
//'   `xmin`, `ymax`, `xmax`, and  `tumoral cells`. `ymin`,
//'   `xmin`, `ymax`, `xmax` report the boundaries of the sampled
//'   rectangular region, while `tumoral cells` is the number of
//'   tumoral cells in the sample.
//' }
//' @field get_species_info Gets the species data\itemize{
//' \item \emph{Returns:} A data frame reporting `mutant` and `epistate`
//'   for each registered species.
//' }
//' @field get_sticks Compute the forest sticks \itemize{
//' \item \emph{Returns:} The list of the forest sticks. Each stick is represented as 
//'   the list of cell identifiers labelling the nodes in the stick
//'   from the higher to the deeper in the forest.
//' }
//' @field get_subforest_for Build a subforest using as leaves some of the original samples \itemize{
//' \item \emph{Parameter:} \code{sample_names} - The names of the samples whose cells will be used
//'   as leaves of the new forest.
//' \item \emph{Returns:} A samples forest built on the samples mentioned in `sample_names`.
//' }
//' @field save Save a phylogenetic forest in a file \itemize{
//' \item \emph{Parameter:} \code{filename} - The path of the file in which the phylogenetic 
//'   forest must be saved.
//' }
  class_<PhylogeneticForest>("PhylogeneticForest")

//' @name PhylogeneticForest$get_nodes
//' @title Get the nodes of the forest
//' @return A data frame representing, for each node
//'   in the forest, the identified (column `cell_id`),
//'   whenever the node is not a root, the ancestor
//'   identifier (column `ancestor`), whenever the
//'   node was sampled, i.e., it is one of the forest
//'   leaves, the name of the sample containing the
//'   node, (column `sample`), the mutant (column
//'   `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' @seealso [SamplesForest$get_nodes()] for usage examples
    .method("get_nodes", (List (PhylogeneticForest::*)() const)(&PhylogeneticForest::get_nodes),
            "Get the nodes of the forest")

//' @name PhylogeneticForest$get_coalescent_cells
//' @title Retrieve most recent common ancestors
//' @description This method retrieves the most recent common ancestors
//'   of a set of cells. If the optional parameter `cell_ids` is
//'   used, this method find the most recent common ancestors of
//'   the cells having an identifier among those in `cell_ids`.
//'   If, otherwise, the optional parameter is not used, this
//'   method find the most recent common ancestors of the forest
//'   leaves.
//' @param cell_ids The list of the identifiers of the cells whose
//'   most recent common ancestors are aimed (optional).
//' @return A data frame representing, for each of the identified
//'   cells, the identified (column `cell_id`), whenever the
//'   node is not a root, the ancestor identifier (column
//'   `ancestor`), whenever the node was sampled, i.e., it is
//'   one of the forest leaves, the name of the sample
//'   containing the node, (column `sample`), the mutant
//'   (column `mutant`), the epistate (column `epistate`),
//'   and the birth time (column `birth_time`).
//' @seealso [SamplesForest$get_coalescent_cells()] for usage examples
    .method("get_coalescent_cells",
            (List (PhylogeneticForest::*)(const std::list<Races::Mutants::CellId>&) const)
                (&PhylogeneticForest::get_coalescent_cells),
            "Get the most recent common ancestor of some cells")
    .method("get_coalescent_cells",
            (List (PhylogeneticForest::*)() const)(&PhylogeneticForest::get_coalescent_cells),
            "Get the most recent common ancestor of all the forest trees")

//' @name PhylogeneticForest$get_subforest_for
//' @title Build a subforest using as leaves some of the original samples
//' @param sample_names The names of the samples whose cells will be used
//'   as leaves of the new forest
//' @return A samples forest built on the samples mentioned in `sample_names`
//' @seealso [SamplesForest$get_subforest_for()] for usage examples
    .method("get_subforest_for", &PhylogeneticForest::get_subforest_for,
            "Get the sub-forest for some of the original samples")

//' @name PhylogeneticForest$get_samples_info
//' @title Retrieve information about the samples
//' @description This method retrieves information about
//'   the samples whose cells were used as leaves
//'   of the samples forest.
//' @return A data frame reporting, for each sample, the
//'   name, the sampling time, the position, and
//'   the number of tumoural cells.
//' @seealso [SamplesForest$get_samples_info()] for usage examples
    .method("get_samples_info", &PhylogeneticForest::get_samples_info,
            "Get some pieces of information about the samples")

//' @name PhylogeneticForest$get_species_info
//' @title Gets the species
//' @return A data frame reporting `mutant` and `epistate`
//'   for each registered species.
    .method("get_species_info", &PhylogeneticForest::get_species_info,
            "Get the recorded species")

//' @name PhylogeneticForest$get_germline_subject
//' @title Gets the germline subject name
//' @return The name of the subject whose germline is used.
    .method("get_germline_subject", &PhylogeneticForest::get_germline_subject,
            "Get the germline subject name")

//' @name PhylogeneticForest$get_sampled_cell_CNAs
//' @title Gets the CNAs of the sampled cells
//' @description This method returns a data frame representing all the CNAs 
//'   in the cells sampled during the simulation and represented by 
//'   the leaves of the phylogenetic forest.
//' @param cell_id The identifier of the cell whose CNAs are aimed (optional).
//' @return A data frame reporting `cell_id`, `type` (`"A"` for amplifications
//'   and `"D"` for deletions), `chromosome`, `begin` (i.e., the first
//'   CNA locus in the chromosome), `end` (i.e., last CNA locus in the
//'   chromosome), `allele`, and `src allele` (the allele origin for
//'   amplifications, `NA` for deletions).
//' @seealso `vignette("mutations")` for usage examples
    .method("get_sampled_cell_CNAs", (List (PhylogeneticForest::*)(const Races::Mutants::CellId&) const)
                (&PhylogeneticForest::get_sampled_cell_CNAs),
            "Get the CNAs of a sampled cell")
    .method("get_sampled_cell_CNAs", (List (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_sampled_cell_CNAs),
            "Get the CNAs of all the sampled cells")

//' @name PhylogeneticForest$get_sampled_cell_SNVs
//' @title Gets the SNVs of the sampled cells
//' @description This method returns a data frame representing all the SNVs 
//'   in the cells sampled during the simulation and represented by 
//'   the leaves of the phylogenetic forest.
//'   The data frame also reports the allele in which SNVs occur to 
//'   support double occurrencies due to CNAs.
//' @param cell_id The identifier of the cell whose SNVs are aimed (optional).
//' @return A data frame reporting `cell_id`, `chromosome`, `chr_pos` (i.e., 
//'   the position in the chromosome), `allele` (in which the SNV
//'   occurs), `ref`, `alt`, and `cause` for each SNV in the sampled
//'   cell genomes.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_sampled_cell_SNVs", (List (PhylogeneticForest::*)(const Races::Mutants::CellId&) const)
                (&PhylogeneticForest::get_sampled_cell_SNVs),
            "Get the SNVs of a sampled cell")
    .method("get_sampled_cell_SNVs", (List (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_sampled_cell_SNVs),
            "Get the SNVs of all the sampled cells")

//' @name PhylogeneticForest$get_germline_SNVs
//' @title Gets the SNVs of the germline
//' @description This method returns a data frame representing all the germinal 
//'   SNVs of the cells represented in the phylogenetic forest.
//'   The data frame also reports the allele in which SNVs occur to 
//'   support double occurrencies due to CNAs.
//' @return A data frame reporting `chromosome`, `chr_pos` (i.e., 
//'   the position in the chromosome), `allele` (in which the SNV
//'   occurs), `ref`, `alt`, and `cause`.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_germline_SNVs", &PhylogeneticForest::get_germline_SNVs,
            "Get the SNVs of the germline")

//' @name PhylogeneticForest$get_sticks
//' @title Compute the forest sticks
//' @description A _crucial node_ of a forest is a root of the forest, a node
//'   whose parent belongs to a different mutant, or the most recent 
//'   common ancestor of two crucial nodes.
//'
//'   A _stick_ is a path of the forest in which the only crucial 
//'   nodes are the first and the last.
//'
//'   This method return the list of the forest sticks. Each stick is
//'   represented by the sequence of cell identifiers labelling the
//'   nodes in the stick.
//' @return The list of the forest sticks. Each stick is represented as 
//'   the list of cell identifiers labelling the nodes in the stick
//'   from the higher to the deeper in the forest.
//' @seealso [SamplesForest$get_sticks()] for usage examples
    .method("get_sticks", (std::list<std::list<Races::Mutants::CellId>> (PhylogeneticForest::*)() const)(&PhylogeneticForest::get_sticks),
            "Get the forest sticks")

//' @name PhylogeneticForest$get_exposures
//' @title Gets the timed exposure data frame
//' @description This method returns a data frame representing the exposure 
//'   evolution over time.
//' @return A data frame reporting `time`, `signature`, `exposure` and, 
//'   `type`.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_exposures", &PhylogeneticForest::get_timed_exposures,
            "Get the timed exposure data frame")

//' @name PhylogeneticForest$get_first_occurrences
//' @title Gets the identifier of the cell in which a mutation occurs for the first time
//' @param mutation A mutation being either a SNV or a CNA.
//' @return The identifier of the cell in which a mutation occurs for the first time.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_first_occurrences", (Rcpp::List (PhylogeneticForest::*)(const SEXP&) const)
                (&PhylogeneticForest::get_first_occurrence),
            "Get the identifier of the cell in which the mutation occurs for the first time")

//' @name PhylogeneticForest$save
//' @title Save a phylogenetic forest in a file
//' @param filename The path of the file in which the phylogenetic 
//'   forest must be saved.
    .method("save", &PhylogeneticForest::save,
            "Save a phylogenetic forest")

    // show
    .method("show", &PhylogeneticForest::show,
            "Describe the PhylogeneticForest");

//' @name load_phylogenetic_forest
//' @title Load a phylogenetic forest from a file
//' @param filename The path of the file from which the phylogenetic 
//'   forest must be load.
//' @return The load phylogenetic forest
  function("load_phylogenetic_forest", &PhylogeneticForest::load,
           "Recover a phylogenetic forest");
}

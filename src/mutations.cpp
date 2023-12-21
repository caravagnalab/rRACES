/*
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023 Alberto Casagrande <alberto.casagrande@uniud.it>
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
//' snv <- SNV("X", 20002, "GAC", "T")
//'
//' # get the chromosome in which `snv` occurs (i.e., "X")
//' snv$get_chromosome()
    .method("get_chromosome",&SNV::get_chromosome, "Get the chromosome of the SNV position")

//' @name SNV$get_position_in_chromosome
//' @title Get the position in chromosome where the SNV occurs.
//' @return The position in chromosome where the SNV occurs.
//' @examples
//' snv <- SNV("X", 20002, "GAC", "T")
//'
//' # get the position in chromosome where `snv` occurs (i.e., 20002)
//' snv$get_position_in_chromosome()
    .method("get_position_in_chromosome",&SNV::get_position_in_chromosome,
            "Get the SNV position in the chromosome")

//' @name SNV$get_context
//' @title Get the context in which the SNV occurs.
//' @return The context in which the SNV occurs.
//' @examples
//' snv <- SNV("X", 20002, "GAC", "T")
//'
//' # get the context in which `snv` occurs (i.e., "GAC")
//' snv$get_context()
    .method("get_context",&SNV::get_context, "Get the SNV context")

//' @name SNV$get_mutated_base
//' @title Get the base after the SNV occurs.
//' @return The base after the SNV occurs.
//' @examples
//' snv <- SNV("X", 20002, "GAC", "T")
//'
//' # get the base after `snv` occurs (i.e., "T")
//' snv$get_mutated_base()
    .method("get_mutated_base",&SNV::get_mutated_base, "Get the new base of the SNV")

//' @name SNV$get_cause
//' @title Get the SNV cause.
//' @description Evey SNV produced by RACES/rRACES is associated to a cause depending
//'        on whether it is part of a genomic characterization of a mutant or it is 
//'        caused by a specific SBS profile. This method return such a cause whenever 
//'        it is available.
//' @return The base after the SNV occurs.
//' @examples
//' # let us build a SNV without specifying any cause for it
//' snv <- SNV("X", 20002, "GAC", "T")
//'
//' # get the cause of `snv` (i.e., "NA")
//' snv$get_cause()
//'
//' # we can also build a SNV, specifying a cause for it
//' snv <- SNV("X", 20002, "GAC", "T", cause = "SBS13")
//'
//' # get the cause of `snv` (i.e., "SBS13")
//' snv$get_cause()
    .method("get_cause",&SNV::get_cause, "Get the cause of the SNV")

//' @name SNV$get_dataframe
//' @title Get a data frame representing the SNV.
//' @description This method returns a data frame that represents the SNV and whose
//'        columns are `chromosome`, `pos_in_chr`, `context`, `mutated_base`, and
//'        `cause`.
//' @examples
//' snv <- SNV("X", 20002, "GAC", "T")
//'
//' snv$get_dataframe()
    .method("get_dataframe",&SNV::get_dataframe, "Get a dataframe representing the SNV")
    .method("show",&SNV::show);

//' @name SNV
//' @title Create a SNV amplification.
//' @param chromosome The name of the chromosome in which the SNV occurs.
//' @param pos_in_chr The position in the chromosome where the SNV occurs.
//' @param context The SNV context.
//' @param mutated_base The base after the mutation.
//' @param cause The cause of the SNV (optional).
//' @examples
//' # create a SNV and do not specify the cause
//' snv <- SNV("X", 20002, "GAC", "T")
//' snv
//'
//' # create a SNV with a cause
//' snv <- SNV("X", 20002, "GAC", "T", cause = "SBS1")
//' snv
  function("SNV", &SNV::build_SNV,
           List::create(_["chromosome"], _["pos_in_chr"], _["context"],
                        _["mutated_base"], _["cause"] = ""),
           "Create a single nucleotide variation (SNV)");

//' @name CNA
//' @title A copy number alteration
  class_<CNA>("CNA")
    .constructor()

//' @name CNA$get_chromosome
//' @title Get the chromosome in which the CNA occurs.
//' @return The chromosome in which the CNA occurs.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
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
//' cna <- Amplification("X", 20002, 100, 1, 0)
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
//'        columns are `chromosome`, `pos_in_chr`, `length`, `mutated_base`, `allele`,
//'        `src_allele`, and `type`.
//' @examples
//' # create an amplification CNA
//' amp_cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' amp_cna$get_dataframe()
//'
//' # create a deletion CNA
//' del_cna <- Deletion("Y", 40020, 200, 0)
//'
//' del_cna$get_dataframe()
    .method("get_dataframe",&CNA::get_dataframe, "Get a dataframe representing the CNA")
    .method("show",&CNA::show);

//' @name CNA
//' @title Create a CNA amplification.
//' @param chromosome The name of the chromosome in which the CNA occurs.
//' @param pos_in_chr The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the CNA occurs.
//' @param src_allele The allele from which the region is amplified (optional:
//'               if not specified the CNA is a deletion).
//' @examples
//' # create an amplification CNA
//' cna <- CNA("X", 20002, 100, 1, 0)
//'
//' cna
  function("CNA", &CNA::build_CNA,
           List::create(_["chromosome"], _["pos_in_chr"], _["len"],
                        _["allele"], _["src_allele"] = -1),
           "Create a copy number alteration (CNA)");

//' @name Amplification
//' @title Create a CNA amplification.
//' @param chromosome The name of the chromosome in which the CNA occurs.
//' @param pos_in_chr The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the amplification is placed.
//' @param src_allele The allele from which the region is amplified.
//' @examples
//' # create an amplification CNA
//' cna <- Amplification("X", 20002, 100, 1, 0)
//'
//' cna
  function("Amplification", &CNA::build_CNA,
           List::create(_["chromosome"], _["pos_in_chr"], _["len"],
                        _["allele"], _["src_allele"]),
           "Create a CNA amplification");

//' @name Deletion
//' @title Create a CNA deletion.
//' @param chromosome The name of the chromosome in which the CNA occurs.
//' @param pos_in_chr The position in the chromosome where the CNA occurs.
//' @param len The CNA length.
//' @param allele The allele in which the deletion occurs.
//' @examples
//' # create a deletion CNA
//' cna <- Deletion("Y", 40020, 200, 0)
//'
//' cna
  function("Deletion", &CNA::build_CNA,
           List::create(_["chromosome"], _["pos_in_chr"], _["len"],
                        _["allele"]),
           "Create a CNA deletion");

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
//'   `MutationEngine$add_coefficients` and `MutationEngine$add_coefficients`
//'   These data are provided by means of the `MutationEngine$add_mutant`.
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

//' @name MutationEngine$add_coefficients
//' @title Add SBS coefficients to the mutation engine. 
//' @description This method adds a set of SBS coefficients to the mutation 
//'        engine. The coefficients must sum up to one.
//'        
//'        The SBS coefficients will be used to establish the probability 
//'        for a SNV to occur depending on its context.
//'
//'        Each set of SBS coefficients is associated to a time that is 
//'        the simulated time in which the set is adopted.
//'        If a time is provided the set of coefficients are used from the
//'        specified time on up to the successive possible change of
//'        coefficients. When a set of SBS coefficients is added to the 
//'        mutation engine without specifying the time, its time is 0.
//' @param time The first simulated time in which the SBS coefficients must 
//'          be taken into account (optional).
//' @param SBS_coefficients A list of SBS coefficients. The coefficients
//'          must sum up to one and their names must corresponds to some of
//'          the SBS name provided in the SBS file.
//' @examples
//' # create a demostrative mutation engine
//' m_engine <- build_mutation_engine("demo")
//'
//' # add a default set of SBS coefficients that will be used from simulated
//' # time 0 up to the successive SBS coefficient change.
//' m_engine$add_coefficients(c(SBS13 = 0.3, SBS1 = 0.7))
//'
//' # add a default set of SBS coefficients that will be used from simulated
//' # time 3.2 up to the end of the simulation.
//' m_engine$add_coefficients(3.2, c(SBS5 = 0.3, SBS2 = 0.2, SBS3 = 0.5))
//'
//' m_engine
    .method("add_coefficients", (void (MutationEngine::*)(const List&))(
              &MutationEngine::add_coefficients), "Add mutational coefficients")
    .method("add_coefficients", (void (MutationEngine::*)(const double&, const List&))(
              &MutationEngine::add_coefficients), "Add mutational coefficients")

//' @name MutationEngine$add_mutant
//' @title Add a mutant specification to the mutation engine.
//' @description This method adds a mutant specification to the mutation engine.
//'           The users must use it to specify the name and the genomic
//'           characterization (i.e., SNVs and CNAs) of all the simulated mutants
//'           together with the mutation rates of its species.
//' @param mutant_name The mutant name.
//' @param species_rates The named list of the species rates whose names are the 
//'           epigenetic states of the species or a single rate, if the mutant
//'           does not have an epigenetic state.
//' @param mutant_SNVs The list of the SNVs characterizing the mutant.
//' @param mutant_CNAs The list of the CNAs characterizing the mutant.
//' @examples
//' # create a demostrative mutation engine
//' m_engine <- build_mutation_engine("demo")
//'
//' # add the mutant "A" characterized by two SNV on chromosome 22 (no CNA) and
//' # having two epigenetic states. Its species "A+" and "A-" have mutation rate 
//' # 3e-9 and 6e-9, respectively.
//' m_engine$add_mutant("A", c("+" = 3e-9, "-" = 6e-9),
//'                     c(SNV("22", 2001, "GAC", "T"),
//'                       SNV("22", 2020, "CCC", "A")))
//'
//' # add the mutant "B" characterized by one SNV on chromosome 1 (no CNA) and
//' # missing epigenetic state. Its species "B" has mutation rate 5e-9.
//' m_engine$add_mutant("B", 5e-9, c(SNV("1", 30, "AGC", "T")))
    .method("add_mutant", (void (MutationEngine::*)(const std::string&, const Rcpp::List& species_rate,
                                                    const Rcpp::List&))(
                                                        &MutationEngine::add_mutant),                                            
            "Add mutant")
    .method("add_mutant", (void (MutationEngine::*)(const std::string&, const Rcpp::List& species_rate,
                                                    const Rcpp::List&, const Rcpp::List&))(
                                                        &MutationEngine::add_mutant), 
            "Add mutant")

//' @name MutationEngine$place_mutations
//' @title Place the mutations on a samples forest
//' @description This method labels each node of a samples forest
//'         by the mutations occuring for the first time in the
//'         cell represented by the node itself and produces a
//'         phylogenetic forest.
//' @param samples_forest A samples forest.
//' @return A phylogenetic forest whose structure corresponds to
//'         `samples_forest`.
//' @examples
//' # create a simulation
//' sim <- new(Simulation)
//' sim$add_mutant("A", 0.2, 0.01)
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
//' m_engine$add_mutant("A", 3e-9, c(SNV("22", 2001, "GAC", "T")))
//'
//' # add the default set of SBS coefficients
//' m_engine$add_coefficients(c(SBS13 = 0.3, SBS1 = 0.7))
//'
//' # place the mutations on the samples forest
//' phylogenetic_forest <- m_engine$place_mutations(samples_forest)
//'
//' phylogenetic_forest
    .method("place_mutations", (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest))(
                                                        &MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")
    .method("place_mutations", (PhylogeneticForest (MutationEngine::*)(const SamplesForest& forest, const int seed))(
                                                        &MutationEngine::place_mutations),
            "Place mutations on a SamplesForest")
    .method("show", &MutationEngine::show);

//' @name build_mutation_engine
//' @title Create a mutation engine object
//' @details This function downloads and sets up the data requires by
//'       by a mutation engine object and build it.
//'       
//'       There are two building modalities: the first one is more general, 
//'       but it requires to specify all the data sources; the second one  
//'       adopts some pre-set configurations, but it is more convient than
//'       the former in many cases.
//'
//'       The first building modality requires to specify the directory in
//'       which the data must be saved, the reference sequence URL, the SBS
//'       file URL, the default number of alleles per chromosome, and, 
//'       optionally, the exceptions on the number of alleles (e.g., one 
//'       allele for "X" and "Y" in human genome) through the parameters 
//'       `directory`, `reference_url`, `SBS_url`, `default_num_of_alleles`,
//'       and `exceptions_on_allele_number`, respectively.
//'
//'       The second building modality exclusively requires a set-up code 
//'       (parameter `setup_code`). The list of supported set-up codes can 
//'       be obtained by using the function `get_mutation_engine_codes`.
//'
//'       The number of context sampling is an optional parameter that allows
//'       sampling the reference contexts while building the context index.
//'       This parameter, which is set to 100 by default, specifies how many
//'       occurences of the same context must be identified before adding 
//'       one of them to the context index. The larger the number of context
//'       sampling, the larger the context index. On the other side, the
//'       lower the number of context sampling, the lower the number of sites
//'       in the reference genome that can be affected by simulated
//'       mutations.
//'
//'       If the parameters of a mutation engine construction match those 
//'       of a previous construction, then the corresponding reference
//'       sequence, the SBS file, and the previously built context index
//'       are loaded from the set-up directory avoiding further 
//'       computations.
//' @seealso `get_mutation_engine_codes()` provides a list of the supported 
//'         set-up codes.
//' @export
//' @param setup_code The set-up code (alternative to `directory`).
//' @param directory The set-up directory (alternative to `setup_code`).
//' @param reference_url The reference genome URL (mandatory when `directory`
//'                is provided).
//' @param SBS_url The SBS file URL (mandatory when `directory` is provided).
//' @param default_num_of_alleles The default number of alleles per chromosome 
//'                (mandatory when `directory` is provided).
//' @param exceptions_on_allele_number A list of exceptions for the default 
//'                number of alleles per chromosome (optional when `directory`
//'                in provided).
//' @param context_sampling The number of reference contexts per context in the
//'     index (optional: default value is 100).
//' @return A mutation engine object.
//' @examples
//' # set the reference and SBS URLs
//' reference_url <- paste0("https://ftp.ensembl.org/pub/release-110/fasta/",
//'                         "homo_sapiens/dna/Homo_sapiens.GRCh38.dna.",
//'                         "chromosome.22.fa.gz")
//' sbs_url <- paste0("https://cancer.sanger.ac.uk/signatures/documents/2124/",
//'                   "COSMIC_v3.4_SBS_GRCh38.txt")
//'
//' # build a mutation engine and save the required files into the
//' # directory "Test"
//' m_engine <- build_mutation_engine(directory = "Test",
//'                                   reference_url = reference_url,
//'                                   SBS_url = sbs_url,
//'                                   default_num_of_alleles = 2)
//'
//' # if the parameters of a mutation engine construction match those of a
//' # previous construction, then the corresponding reference sequence,
//' # the SBS file, and the previously built context index are loaded from
//' # the set-up directory avoiding further computations.
//' m_engine <- build_mutation_engine("Test", reference_url, sbs_url,
//'                                   default_num_of_alleles = 2)
//'
//' # if the `context_sampling` parameter changes, a new context index is
//' # built, while neither the reference sequence nor the SBS file are 
//' # downloaded again.
//' m_engine <- build_mutation_engine("Test", reference_url, sbs_url,
//'                                   default_num_of_alleles = 2,
//'                                   context_sampling = 50)
//'
//' # a futher contruction with the same parameters avoids both downloads
//' # and context index construction.
//' m_engine <- build_mutation_engine("Test", reference_url, sbs_url,
//'                                   default_num_of_alleles = 2,
//'                                   context_sampling = 50)
//'
//' # the parameters `directory`, `reference_url`, and `SBS_url` can be
//' # avoided by providing the `setup_code` parameter. The set-up code 
//' # `demo` is provided among those available for testing purpose.
//' m_engine <- build_mutation_engine(setup_code = "demo")
//' 
//' # the `context_sampling` can be used also when a pre-defined set-up
//' # configuration is adopted.
//' m_engine <- build_mutation_engine(setup_code = "demo",
//'                                   context_sampling = 50)
//' 
//' # remove the "Test" and "demo" directories
//' unlink("Test", recursive = TRUE)
//' unlink("demo", recursive = TRUE)
  function("build_mutation_engine", &MutationEngine::build_MutationEngine,
           List::create(_["directory"] = "",
                        _["reference_url"] = "",  _["SBS_url"] = "",
                        _["default_num_of_alleles"] = 0,
                        _["exceptions_on_allele_number"] = Rcpp::List::create(),
                        _["setup_code"] = "", 
                        _["context_sampling"] = 100),
           "Create a MutationEngine");

//' @name get_mutation_engine_codes
//' @title Get the supported codes for predefined set-up
//' @return A data frame reporting the code and a description for each
//'      supported predefined set-up.
//' @seealso `build_mutation_engine()` to build a mutation engine
//' @export
//' @examples
//' # get the list of supported mutation engine set-up codes
//' get_mutation_engine_codes()
  function("get_mutation_engine_codes", &MutationEngine::get_supported_setups,
           "Get mutation engine supported codes");

//' @name PhylogeneticForest
//' @title The phylogenetic forest of cells in samples.
//' @description Represents the phylogenetic forest of the
//'       cells sampled during the computation. The leaves of
//'       this forest are the sampled cells.
//'       This class is analoguous to the class `SamplesForest`, 
//'       but each node is labelled with the mutations occuring
//'       for the first time on the cell represented by the node
//'       itself. Moreover each leaf is also associated with the
//'       genome mutations occurring in the corresponding cell.
//' @field get_coalescent_cells Retrieve most recent common ancestors\itemize{
//' \item \emph{Parameter:} \code{cell_ids} - The list of the identifiers of the
//'               cells whose most recent common ancestors are aimed (optional).
//' \item \emph{Return:} A data frame representing, for each of the identified
//'         cells, the identified (column "cell_id"), whenever the
//'         node is not a root, the ancestor identifier (column
//'         "ancestor"), whenever the node was sampled, i.e., it is
//'         one of the forest leaves, the name of the sample
//'         containing the node, (column "sample"), the mutant
//'         (column "mutant"), the epistate (column "epistate"),
//'         and the birth time (column "birth_time").
//' }
//' @field get_nodes Get the forest nodes \itemize{
//' \item \emph{Return:} A data frame representing, for each node
//'              in the forest, the identified (column "id"),
//'              whenever the node is not a root, the ancestor
//'              identifier (column "ancestor"), whenever the node
//'              was sampled, i.e., it is one of the forest
//'              leaves, the name of the sample containing the
//'              node, (column "sample"), the mutant (column
//'              "mutant"), the epistate (column "epistate"),
//'              and the birth time (column "birth_time").
//' }
//' @field get_samples_info Retrieve information about the samples \itemize{
//' \item \emph{Returns:} A data frame containing, for each sample collected
//'         during the simulation, the columns "name", "time", "ymin",
//'         "xmin", "ymax", "xmax", and  "tumoral cells". "ymin",
//'         "xmin", "ymax", "xmax" report the boundaries of the sampled
//'         rectangular region, while "tumoral cells" is the number of
//'         tumoral cells in the sample.
//' }
//' @field get_sampled_cell_SNVs Gets a the SNVs of the sampled cells \itemize{
//' \item \emph{Returns:} A data frame reporting `cell_id`, `chromosome`, 
//'          `chr_pos` (i.e., position in the chromosome), `allele` (in which
//'          the SNV occurs), `context`, `mutated_base`, and `cause` for each
//'          SNV in the sampled cell genomes.
//' }
//' @field get_species_info Gets the species data\itemize{
//' \item \emph{Returns:} A data frame reporting "mutant" and "epistate"
//'            for each registered species.
//' }
//' @field get_subforest_for Build a subforest using as leaves some of the original samples \itemize{
//' \item \emph{Parameter:} \code{sample_names} - The names of the samples whose cells will be used
//'         as leaves of the new forest.
//' \item \emph{Returns:} A samples forest built on the samples mentioned in `sample_names`.
//' }
//' @field save Save a phylogenetic forest in a file \itemize{
//' \item \emph{Parameter:} \code{filename} - The path of the file in which the phylogenetic 
//'            forest must be saved.
//' }
  class_<PhylogeneticForest>("PhylogeneticForest")

//' @name PhylogeneticForest$get_nodes
//' @title Get the nodes of the forest
//' @return A data frame representing, for each node
//'         in the forest, the identified (column "cell_id"),
//'         whenever the node is not a root, the ancestor
//'         identifier (column "ancestor"), whenever the
//'         node was sampled, i.e., it is one of the forest
//'         leaves, the name of the sample containing the
//'         node, (column "sample"), the mutant (column
//'         "mutant"), the epistate (column "epistate"),
//'         and the birth time (column "birth_time").
//' @seealso `SamplesForest$get_nodes()` for usage examples
    .method("get_nodes", (List (PhylogeneticForest::*)() const)(&PhylogeneticForest::get_nodes),
            "Get the nodes of the forest")

//' @name PhylogeneticForest$get_coalescent_cells
//' @title Retrieve most recent common ancestors
//' @description This method retrieves the most recent common ancestors
//'         of a set of cells. If the optional parameter `cell_ids` is
//'         used, this method find the most recent common ancestors of
//'         the cells having an identifier among those in `cell_ids`.
//'         If, otherwise, the optional parameter is not used, this
//'         method find the most recent common ancestors of the forest
//'         leaves.
//' @param cell_ids The list of the identifiers of the cells whose
//'         most recent common ancestors are aimed (optional).
//' @return A data frame representing, for each of the identified
//'         cells, the identified (column "cell_id"), whenever the
//'         node is not a root, the ancestor identifier (column
//'         "ancestor"), whenever the node was sampled, i.e., it is
//'         one of the forest leaves, the name of the sample
//'         containing the node, (column "sample"), the mutant
//'         (column "mutant"), the epistate (column "epistate"),
//'         and the birth time (column "birth_time").
//' @seealso `SamplesForest$get_coalescent_cells()` for usage examples
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
//'         as leaves of the new forest
//' @return A samples forest built on the samples mentioned in `sample_names`
//' @seealso `SamplesForest$get_subforest_for()` for usage examples
    .method("get_subforest_for", &PhylogeneticForest::get_subforest_for,
            "Get the sub-forest for some of the original samples")

//' @name PhylogeneticForest$get_samples_info
//' @title Retrieve information about the samples
//' @description This method retrieves information about
//'           the samples whose cells were used as leaves
//'           of the samples forest.
//' @return A data frame reporting, for each sample, the
//'           name, the sampling time, the position, and
//'           the number of tumoural cells.
//' @seealso `SamplesForest$get_samples_info()` for usage examples
    .method("get_samples_info", &PhylogeneticForest::get_samples_info,
            "Get some pieces of information about the samples")

//' @name PhylogeneticForest$get_species_info
//' @title Gets the species
//' @return A data frame reporting "mutant" and "epistate"
//'            for each registered species.
    .method("get_species_info", &PhylogeneticForest::get_species_info,
            "Get the recorded species")

//' @name PhylogeneticForest$get_sampled_cell_SNVs
//' @title Gets a the SNVs of the sampled cells
//' @description This method returns a data frame representing all the SNVs 
//'          in the cells sampled during the simulation and represented by 
//'          the leaves of the phylogenetic forest.
//'          The data frame also reports the allele in which SNVs occur to 
//'          support double occurrencies due to CNAs.
//' @param cell_id The identifier of the cell whose SNVs are aimed (optional).
//' @return A data frame reporting `cell_id`, `chromosome`, `chr_pos` (i.e., 
//'          position in the chromosome), `allele` (in which the SNV occurs),
//'          `context`, `mutated_base`, and `cause` for each SNV in the 
//'          sampled cell genomes.
//' @seealso `vignette("mutations")` for usage examples
    .method("get_sampled_cell_SNVs", (List (PhylogeneticForest::*)(const Races::Mutants::CellId&) const)
                (&PhylogeneticForest::get_sampled_cell_SNVs),
            "Get the SNVs of a sampled cell")
    .method("get_sampled_cell_SNVs", (List (PhylogeneticForest::*)() const)
                (&PhylogeneticForest::get_sampled_cell_SNVs),
            "Get the SNVs of all the sampled cells")

//' @name PhylogeneticForest$save
//' @title Save a phylogenetic forest in a file
//' @param filename The path of the file in which the phylogenetic 
//'            forest must be saved.
    .method("save", &PhylogeneticForest::save,
            "Save a phylogenetic forest")

    // show
    .method("show", &PhylogeneticForest::show,
            "Describe the PhylogeneticForest");

//' @name load_phylogenetic_forest
//' @title Load a phylogenetic forest from a file
//' @param filename The path of the file from which the phylogenetic 
//'            forest must be load.
//' @return The load phylogenetic forest
  function("load_phylogenetic_forest", &PhylogeneticForest::load,
           "Recover a phylogenetic forest");
}

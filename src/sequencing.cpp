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

#include "sequencers.hpp"
#include "seq_simulation.hpp"
#include "sampled_cell.hpp"

using namespace Rcpp;

RCPP_MODULE(Sequencing){

//' @name ErrorlessIlluminaSequencer
//' @title An error-less Illumina sequencer class
//' @description This class implements a perferct Illumina sequencers that
//'   does not commit errors.
//' @seealso `simulate_seq()`, `simulate_normal_seq()`, and
//'   `vignette("sequencing")` for usage examples
  class_<ErrorlessIlluminaSequencer>("ErrorlessIlluminaSequencer")

    .method("show", &ErrorlessIlluminaSequencer::show,
            "Show a description for the sequencer");

//' @name ErrorlessIlluminaSequencer
//' @description This method builds an error-less Illumina 
//'   sequencer model.
//' @return A new error-less Illumina sequencer.
//' @examples
//' # build a sequencer model
//' sequencer <- ErrorlessIlluminaSequencer()
//' sequencer
    function("ErrorlessIlluminaSequencer",
             &ErrorlessIlluminaSequencer::build_sequencer,
             "Build a new error-less Illumina sequencer");

//' @name BasicIlluminaSequencer
//' @title A basic Illumina sequencer class
//' @description This class implements a basic model for Illumina sequencers.
//' @details It specifies a simulated sequencing error rate and the simulated
//'   sequencing errors will occurs according to that rate.
//' @seealso `simulate_seq()`, `simulate_normal_seq()`, and
//'   `vignette("sequencing")` for usage examples
  class_<BasicIlluminaSequencer>("BasicIlluminaSequencer")
    .method("show", &BasicIlluminaSequencer::show,
            "Show a description for the sequencer")

//' @name BasicIlluminaSequencer$get_error_rate
//' @title Getting error rate
//' @description This method returns the sequencing error rate of the
//'   simulated illumina sequencer.
//' @return The sequencing error rate of the simualted sequencer.
//' @examples
//' # build a basic Illumina sequencer model whose errors occur
//' # at rate 4e-3
//' sequencer <- BasicIlluminaSequencer(4e-3)
//'
//' sequencer$get_error_rate()
    .method("get_error_rate",
            (const double& (BasicIlluminaSequencer::*)())(
                &BasicIlluminaSequencer::get_error_rate),
            "Get the sequencer error rate");

//' @name BasicIlluminaSequencer
//' @description This method builds a basic Illumina sequencer model.
//' @param error_rate The error rate of the sequencer model.
//' @param random_quality_scores A Boolean flag to enable a basic quality score 
//'   model. When it is set to `false`, all the bases with no sequencing
//'   errors have the same quality score. The random quality score model
//'   increases the computation time of about 70%. (default: `true`)
//' @param seed The seed of the random error generator (default: `NULL`).
//' @return A basic Illumina sequencer model.
//' @examples
//' # build a sequencer model having error rate 4e-3
//' sequencer <- BasicIlluminaSequencer(error_rate=4e-3)
//' sequencer
//'
//' # build a sequencer model having error rate 4e-3 and set the seed to 5
//' sequencer <- BasicIlluminaSequencer(error_rate=4e-3, seed=5)
//' sequencer
    function("BasicIlluminaSequencer", &BasicIlluminaSequencer::build_sequencer,
             List::create(_["error_rate"], _["random_quality_scores"] = true,
                          _["seed"] = R_NilValue),
             "Create a basic Illumina sequencer model");

//' @name simulate_seq
//' @title Simulating the sequencing
//' @description This method simulates the sequencing of the samples in a phylogenetic
//'   forest.
//' @param phylo_forest A phylogenetic forest.
//' @param sequencer The sequencer that performs the sequencing simulation
//'   (default: an `ErrorlessIlluminaSequencer`).
//' @param chromosomes The chromosomes that must be considered (default:
//'   `NULL`, i.e., all the reference chromosomes).
//' @param coverage The sequencing coverage (default: `10`).
//' @param read_size The read size (default: `150`).
//' @param insert_size_mean The insert size mean. Use 0 for single read 
//'   sequencing and any value greater than 0 for pair read sequencing
//'   (default: `0`).
//' @param insert_size_stddev The insert size standard deviation.
//'   (default: `10`).
//' @param output_dir The SAM output directory (default:
//'   `"rRACES_SAM"`).
//' @param write_SAM A Boolean flag to enable/disable SAM generation
//'   (default: `FALSE`).
//' @param update_SAM Update the output directory (default: `FALSE`).
//' @param cell_labelling The labelling function for sampled cells
//'   See `vignette("sample_partition")` for details (default: `NULL`).
//' @param purity The ratio between the number of sample tumeral cell
//'   and that of all the cells, i.e., tumour and normal
//'   ones. This value must belong to the interval [0,1]
//'   (default: `1`).
//' @param with_normal_sample A Boolean flag to enable/disable the
//'   analysis of a normal sample (default: `TRUE`).
//' @param filename_prefix The prefix of the output SAM file name
//'   (default: `"chr_"`).
//' @param template_name_prefix The template name prefix (default:
//'   `"r"`).
//' @param seed The random seed for the internal random generator
//'   (optional).
//' @return A dataframe representing, for each of the observed SNVs
//'   and indels, the chromosome and the position in which it occurs
//'   (columns `chr` and `chr_pos`), the mutation reference and
//'   alterate sequences (columns `ref` and `alt`, respectively), its
//'   cause and class (columns `causes`, and `classes`, respectively).
//'   Moreover, for each of the sequencied samples `<sample name>`,
//'   the returned dataframe contains three columns: the number of
//'   reads in which the corresponding mutation occurs (column
//'   `<sample name>.occurrences`), the coverage of the mutation
//'   (column `<sample name>.coverage`), and the corresponding VAF
//'   (column `<sample name>.VAF`).
//' @seealso `BasicIlluminaSequencer` and
//'   `ErrorlessIlluminaSequencer` as sequencer types, and
//'   `vignette("sequencing")` for usage examples
  function("simulate_seq", &simulate_seq,
           List::create(_["phylo_forest"], _["sequencer"] = R_NilValue,
                        _["chromosomes"] = R_NilValue,
                        _["coverage"] = 10,
                        _["read_size"] = 150, _["insert_size_mean"] = 0,
                        _["insert_size_stddev"] = 10, 
                        _["output_dir"] = "rRACES_SAM",
                        _["write_SAM"] = false, _["update_SAM"] = false,
                        _["cell_labelling"] = R_NilValue, _["purity"] = 1,
                        _["with_normal_sample"] = true,
                        _["filename_prefix"] = "chr_",
                        _["template_name_prefix"] = "r",
                        _["seed"] = R_NilValue),
           "Simulate the sequencing of the samples in a phylogenetic forest");

//' @name simulate_normal_seq
//' @title Simulating wild-type sequencing
//' @description This method simulates a wild-type sample sequencing in a phylogenetic
//'   forest. Add the cells in the wild-type sample contains the germline mutations.
//'   The forest pre-neoplastic mutations are also added to the sample by default. 
//'   However, they can be avoided by using the parameter `with_preneoplastic`.
//' @param phylo_forest A phylogenetic forest.
//' @param sequencer The sequencer that performs the sequencing simulation
//'   (default: an `ErrorlessIlluminaSequencer`).
//' @param chromosomes The chromosomes that must be considered (default:
//'   `NULL`, i.e., all the reference chromosomes).
//' @param coverage The sequencing coverage (default: `10`).
//' @param read_size The read size (default: `150`).
//' @param insert_size_mean The insert size mean. Use 0 for single read 
//'   sequencing and any value greater than 0 for pair read sequencing
//'   (default: `0`).
//' @param insert_size_stddev The insert size standard deviation.
//'   (default: `10`).
//' @param output_dir The SAM output directory (default:
//'   `"rRACES_normal_SAM"`).
//' @param write_SAM A Boolean flag to enable/disable SAM generation
//'   (default: `TRUE`).
//' @param update_SAM Update the output directory (default: `FALSE`).
//' @param with_preneoplastic Add the forest pre-neoplastic mutations
//'   to the sample cells. (default: `TRUE`).
//' @param filename_prefix The prefix of the output SAM file name
//'   (default: `"chr_"`).
//' @param template_name_prefix The template name prefix (default:
//'   `"r"`).
//' @param seed The random seed for the internal random generator
//'   (optional).
//' @return A dataframe representing, for each of the observed
//'   SNVs, the chromosome and the position in which
//'   it occurs (columns `chr` and `chr_pos`),
//'   the SNV reference base, the alterate base, the causes,
//'   and the classes of the SNV (columns `ref_base`, `alt_base`,
//'   `causes`, and `classes`, respectively). Moreover, for each
//'   of the sequencied samples `normal_sample`, the returned
//'   dataframe contains three columns: the number of reads in
//'   which the corresponding SNV occurs (column
//'   `normal_sample.occurrences`), the coverage of the SNV
//'   locus (column `normal_sample.coverage`), and the
//'   corresponding VAF (column `normal_sample.VAF`).
//' @seealso `BasicIlluminaSequencer` and
//'   `ErrorlessIlluminaSequencer` as sequencer types, and
//'   `vignette("sequencing")` for usage examples
  function("simulate_normal_seq", &simulate_normal_seq,
           List::create(_["phylo_forest"], _["sequencer"] = R_NilValue,
                        _["chromosomes"] = R_NilValue,
                        _["coverage"] = 10,
                        _["read_size"] = 150, _["insert_size_mean"] = 0,
                        _["insert_size_stddev"] = 10, 
                        _["output_dir"] = "rRACES_normal_SAM",
                        _["write_SAM"] = true, _["update_SAM"] = false,
                        _["with_preneoplastic"] = true,
                        _["filename_prefix"] = "chr_",
                        _["template_name_prefix"] = "r",
                        _["seed"] = R_NilValue),
           "Simulate the sequencing of a normal sample");

//' @name SampledCell
//' @title A sampled cell
//' @description The sampled cell class for sample labelling.
//' @details There is no public constructor for this class as it is
//'   exclusively used by `simulate_seq()` to label sampled cells.
//' @seealso `simulate_seq()` and `vignette("sample_partition")`
  class_<SampledCell>("SampledCell")

//' @name SampledCell$epistate
//' @title Getting the sampled cell epigenetic state
//' @description The epigenetic state of the sampled cell.
//' @details This property is the epigenetic state of the sampled cell.
//'   It can be one among "`+`", "`-`", or "".
//' @seealso `simulate_seq()` and `vignette("sample_partition")`
    .property("epistate", &SampledCell::epistate,
              "The cell epistate")

//' @name SampledCell$mutant
//' @title Getting the sampled cell mutant
//' @description The mutant name of the sampled cell.
//' @details This property is the mutant name of the sampled cell.
//' @seealso `simulate_seq()` and `vignette("sample_partition")`
    .property("mutant", &SampledCell::mutant,
              "The cell mutant name")

//' @name SampledCell$species
//' @title Getting the sampled cell species
//' @description The species name of the sampled cell.
//' @details This property is the species name of the sampled cell.
//' @seealso `simulate_seq()` and `vignette("sample_partition")`
    .property("species", &SampledCell::species,
              "The cell species name")

//' @name SampledCell$birth_time
//' @title Getting the sampled cell birth time
//' @description The birth time of the sampled cell.
//' @details This property is the birth time of the sampled cell.
//' @seealso `simulate_seq()` and `vignette("sample_partition")`
    .property("birth_time", &SampledCell::birth_time,
              "The cell birth time")

//' @name SampledCell$mutations
//' @title Getting the sampled cell mutations
//' @description The mutations of the sampled cell.
//' @details This property contains a dataframe that represents the sampled
//'   cell mutations. The dataframe format is analoguous to that returned by
//'   `PhylogeneticForest$get_sampled_cell_mutations()`: it has columns
//'   `cell_id`, `chr`, (i.e., the mutation chromosome), `chr_pos` (i.e.,
//'   position in the chromosome), `allele` (in which the mutation occurs),
//'   `ref`, `alt`, `type` (i.e., either `"SNV"` or `"indel"`), `cause`, and
//'   `class` (i.e., `"driver"`, `"passenger"`, `"germinal"` or
//'   `"preneoplastic"`).
//' @seealso `simulate_seq()` and `vignette("sample_partition")`
    .property("mutations", &SampledCell::mutations,
              "The cell mutation dataframe");
}

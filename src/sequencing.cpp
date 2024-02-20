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

#include "seq_simulation.hpp"

using namespace Rcpp;

RCPP_MODULE(Sequencing){

//' @name simulate_seq
//' @title Simulate the sequencing of the samples in a phylogenetic forest
//' @param phylo_forest A phylogenetic forest.
//' @param coverage The sequencing coverage (default value: 1).
//' @param read_size The read size (default value: 150).
//' @param insert_size The insert size. Use 0 for single read sequencing 
//'              and any value greater than 0 for pair read sequencing 
//'              (default value: 0).
//' @param output_dir The SAM output directory (default value: 
//'              "rRACES_SAM").
//' @param write_SAM A Boolean flag to enable/disable SAM generation 
//'              (default: FALSE).
//' @param epi_FACS Perform an epigenetic FACS analysis (default: FALSE).
//' @param rnd_seed The random seed for the internal random generator 
//'              (default value: 0).
//' @return A data frame representing, for each of the observed
//'         SNVs, the chromosome and the position in which
//'         it occurs (columns `chromosome` and `chr_pos`),
//'         the SNV context and the new base (columns `context`
//'         and `alt_base`). Moreover, for each of the 
//'         sequencied samples `<sample name>`, the returned data 
//'         frame contains three columns: the number of reads 
//'         in which the corresponding SNV occurs 
//'         (column `<sample name>.occurrences`), the coverage of
//'         the SNV locus (column `<sample name>.coverage`), and the
//'         corresponding VAF (column `<sample name>.VAF`).
//' @seealso `vignette("sequencing")` for usage examples
  function("simulate_seq", &simulate_seq,
           List::create(_["phylo_forest"], _["coverage"]=1,
                        _["read_size"] = 150, _["insert_size"] = 0,
                        _["output_dir"] = "rRACES_SAM",
                        _["write_SAM"] = false, _["epi_FACS"] = false, 
                        _["rnd_seed"] = 0),
           "Simulate the sequencing of the samples in a phylogenetic forest");
}

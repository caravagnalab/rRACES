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

#include "setup_mutation_engine.hpp"

using namespace Rcpp;

RCPP_MODULE(Mutations){
//' @name setup_mutation_engine
//' @title Download and set up mutation engine data
//' @details This function downloads and sets up the data requires by 
//'       by the mutation engine.
//'       The number of context sampling is an optional parameter that allows
//'       sampling the reference contexts while building the context index.
//'       This parameter, which is set to 100 by default, specifies how many
//'        occurences of the same context must be identified before adding one of 
//'       them to the context index. The larger the number of context sampling, 
//'       the larger the context index. On the other side, the lower the number
//'       of context sampling, the lower the number of sites in the reference
//'       genome that can be affected by simulated mutations.
//' @export
//' @param directory The setup directory.
//' @param reference_url The reference genome URL.
//' @param SBS_url The SBS file URL.
//' @param context_sampling The number of reference contexts per context in the 
//'     index (optional: default value is 100).
//' @seealso [set_mutation_engine_by_code()] for an easier-to-use version of this
//'            this function.
//' @examples
//' # set the reference and SBS URLs
//' reference_URL <- paste0("https://ftp.ensembl.org/pub/release-110/fasta/",
//'                         "homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz")
//' SBS_URL <- paste0("https://cancer.sanger.ac.uk/signatures/documents/2124/",
//'                   "COSMIC_v3.4_SBS_GRCh38.txt")
//'
//' # build the set-up and store it in the directory "Test"
//' setup_mutation_engine("Test", reference_URL, SBS_URL)
//'
//' # remove the "Test" directory
//' unlink("Test", recursive=TRUE)
//'
//' # re-build the set-up with `context_sampling` equals to 50
//' setup_mutation_engine("Test", reference_URL, SBS_URL, 
//'                       context_sampling = 50)
//'
//' # remove the "Test" directory
//' unlink("Test", recursive=TRUE)
  function("setup_mutation_engine", (void (*)(const std::string&, const std::string&, 
                                              const std::string&, 
                                              const size_t&))&setup_mutation_engine,
           List::create(_["directory"], _["reference_url"], _["SBS_url"],
                        _["context_sampling"] = 100),
           "Setup mutation engine");

//' @name get_mutation_engine_codes
//' @title Get the supported codes for predefined set-up
//' @return A data frame reporting the code and a description for each 
//'      supported predefined set-up.
//' @seealso [set_mutation_engine_by_code()]
//' @export
//' @examples
//' get_mutation_engine_codes()
  function("get_mutation_engine_codes", &get_mutation_engine_supported_codes,
           "Get mutation engine supported codes");

//' @name setup_mutation_engine_by_code
//' @title Download and set up mutation engine data
//' @details This function downloads and sets up the data requires by 
//'       by the mutation engine.
//'       The number of context sampling is an optional parameter that allows
//'       sampling the reference contexts while building the context index.
//'       This parameter, which is set to 100 by default, specifies how many
//'       occurences of the same context must be identified before adding one of 
//'       them to the context index. The larger the number of context sampling, 
//'       the larger the context index. On the other side, the lower the number
//'       of context sampling, the lower the number of sites in the reference
//'       genome that can be affected by simulated mutations.
//' @export
//' @param setup_code The code of the set-up.
//' @param context_sampling The number of reference contexts per context in the 
//'     index (optional: default value is 100).
//' @seealso [set_mutation_engine()] for a more flexible vesion of this function.
//' @examples
//' # build the predefined set-up "demo"
//' setup_mutation_engine_by_code("demo")
//'
//' # remove the "demo" directory
//' unlink("demo", recursive=TRUE)
//'
//' # build the predefined set-up "demo"
//' setup_mutation_engine_by_code("demo", 50)
//'
//' # remove the "demo" directory
//' unlink("demo", recursive=TRUE)
  function("setup_mutation_engine_by_code", 
           (void (*)(const std::string&, const size_t&))&setup_mutation_engine_by_code,
           List::create(_["setup_code"], _["context_sampling"] = 100),
           "Setup mutation engine with code");
}
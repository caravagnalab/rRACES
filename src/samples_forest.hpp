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

#ifndef __RRACES_SAMPLES_FOREST__
#define __RRACES_SAMPLES_FOREST__

#include <vector>

#include <Rcpp.h>

#include <phylogenetic_forest.hpp>


//' @name SamplesForest
//' @title The forest of the sampled cell ancestors.
//' @description Represents the forest of the ancestors of the
//'       cells sampled during the computation. The leaves of
//'       this forest are the sampled cells.
//' @field get_coalescent_cells Retrieve most recent common ancestors\itemize{
//' \item \emph{Parameter:} \code{cell_ids} - The list of the identifiers of the
//'               cells whose most recent common ancestors are aimed (optional).
//' \item \emph{Return:} A data frame representing, for each of the identified
//'         cells, the identified (column "cell_id"), whenever the
//'         node is not a root, the ancestor identifier (column
//'         "ancestor"), whenever the node was sampled, i.e., it is
//'         one of the forest leaves, the name of the sample
//'         containing the node, (column "sample"), the genotype
//'         (column "genotype"), the epistate (column "epistate"),
//'         and the birth time (column "birth_time").
//' }
//' @field get_nodes Get the forest nodes \itemize{
//' \item \emph{Return:} A data frame representing, for each node
//'              in the forest, the identified (column "id"),
//'              whenever the node is not a root, the ancestor
//'              identifier (column "ancestor"), whenever the node
//'              was sampled, i.e., it is one of the forest
//'              leaves, the name of the sample containing the
//'              node, (column "sample"), the genotype (column
//'              "genotype"), the epistate (column "epistate"),
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
//' @field get_species_info Gets the species data\itemize{
//' \item \emph{Returns:} A data frame reporting "genotype" and "epistate"
//'            for each registered species.
//' }
//' @field get_subforest_for Build a subforest using as leaves some of the original samples \itemize{
//' \item \emph{Parameter:} \code{sample_names} - The names of the samples whose cells will be used
//'         as leaves of the new forest.
//' \item \emph{Returns:} A samples forest built on the samples mentioned in `sample_names`.
//' }
class SamplesForest : private Races::Drivers::DescendantsForest
{
  SamplesForest();

  Rcpp::List get_nodes(const std::vector<Races::Drivers::CellId>& cell_ids) const;

public:
  SamplesForest(const Races::Drivers::Simulation::Simulation& simulation);

  Rcpp::List get_nodes() const;

  Rcpp::List get_samples_info() const;

  Rcpp::List get_species_info() const;

  Rcpp::List get_coalescent_cells() const;

  Rcpp::List get_coalescent_cells(const std::list<Races::Drivers::CellId>& cell_ids) const;

  SamplesForest get_subforest_for(const std::vector<std::string>& sample_names) const;

  void show() const;
};

//' @name SamplesForest$get_nodes
//' @title Get the nodes of the forest
//' @return A data frame representing, for each node
//'         in the forest, the identified (column "cell_id"),
//'         whenever the node is not a root, the ancestor
//'         identifier (column "ancestor"), whenever the
//'         node was sampled, i.e., it is one of the forest
//'         leaves, the name of the sample containing the
//'         node, (column "sample"), the genotype (column
//'         "genotype"), the epistate (column "epistate"),
//'         and the birth time (column "birth_time").
//' @examples
//' # create a simulation having name "get_nodes_test"
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  growth_rate = 0.2,
//'                  death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' forest$get_nodes()


//' @name SamplesForest$get_samples_info
//' @title Retrieve information about the samples
//' @description This method retrieves information about
//'           the samples whose cells were used as leaves
//'           of the samples forest.
//' @return A data frame reporting, for each sample, the
//'           name, the sampling time, the position, and
//'           the number of tumoural cells.
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  growth_rate = 0.2,
//'                  death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' # get information about the sampled whose cells
//' # are the forest leaves, i.e, S1 and S2
//' forest$get_samples_info()


//' @name SamplesForest$get_coalescent_cells
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
//'         containing the node, (column "sample"), the genotype
//'         (column "genotype"), the epistate (column "epistate"),
//'         and the birth time (column "birth_time").
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  growth_rate = 0.2,
//'                  death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' forest$get_coalescent_cells()


//' @name SamplesForest$get_subforest_for
//' @title Build a subforest using as leaves some of the original samples
//' @param sample_names The names of the samples whose cells will be used
//'         as leaves of the new forest
//' @return A samples forest built on the samples mentioned in `sample_names`
//' @examples
//' sim <- new(Simulation)
//' sim$add_genotype(genotype = "A",
//'                  growth_rate = 0.2,
//'                  death_rate = 0.01)
//' sim$place_cell("A", 500, 500)
//'
//' sim$death_activation_level <- 100
//' sim$run_up_to_size(species = "A", num_of_cells = 50000)
//'
//' # sample the region [450,500]x[475,550]
//' sim$sample_cells("S1", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' sim$run_up_to_size(species = "A", num_of_cells = 60000)
//'
//' # sample again the same region
//' sim$sample_cells("S2", lower_corner=c(450,475), upper_corner=c(500,550))
//'
//' # build the samples forest
//' forest <- sim$get_samples_forest()
//'
//' forest$get_subforest_for("S2")


//' @name SamplesForest$get_species_info
//' @title Gets the species
//' @return A data frame reporting "genotype" and "epistate"
//'            for each registered species.


#endif // __RRACES_SAMPLES_FOREST__

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

#include "phylogenetic_forest.hpp"

PhylogeneticForest::PhylogeneticForest():
  Races::Mutations::PhylogeneticForest()
{}

PhylogeneticForest::PhylogeneticForest(const Races::Mutations::PhylogeneticForest& orig):
   Races::Mutations::PhylogeneticForest(orig)
{}

PhylogeneticForest PhylogeneticForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
  PhylogeneticForest forest;

  static_cast<Races::Mutations::PhylogeneticForest&>(forest) = Races::Mutations::PhylogeneticForest::get_subforest_for(sample_names);

  return forest;
}

size_t count_SNVs(const Races::Mutations::CellGenomeMutations& cell_mutations)
{
  size_t counter{0};
  for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {
    for (const auto& [allele_id, allele]: chromosome.get_alleles()) {
      for (const auto& [fragment_pos, fragment]: allele.get_fragments()) {
          counter += fragment.get_SNVs().size();
      }
    }
  }

  return counter;
}

size_t count_SNVs(const std::map<Races::Mutants::CellId, Races::Mutations::CellGenomeMutations>& genome_mutations)
{
  size_t counter{0};
  for (const auto& [cell_id, cell_mutations]: genome_mutations) {
    counter += count_SNVs(cell_mutations);
  }

  return counter;
}

void fill_lists(const Races::Mutations::CellGenomeMutations& cell_mutations,
                Rcpp::IntegerVector& cell_ids, Rcpp::CharacterVector& chr_names, 
                Rcpp::IntegerVector& chr_pos, Rcpp::IntegerVector& alleles,
                Rcpp::CharacterVector& contexts, Rcpp::CharacterVector& mutated_bases,
                Rcpp::CharacterVector& causes, size_t& index)
{
  using namespace Races::Mutations;

  for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {
    for (const auto& [allele_id, allele]: chromosome.get_alleles()) {
      for (const auto& [fragment_pos, fragment]: allele.get_fragments()) {
        for (const auto& [snv_pos, snv]: fragment.get_SNVs()) {
          cell_ids[index] = cell_mutations.get_id();
          chr_names[index] = GenomicPosition::chrtos(chr_id);
          chr_pos[index] = snv.position;
          alleles[index] = allele_id;
          contexts[index] = snv.context.get_sequence();
          mutated_bases[index] = std::string(1,snv.mutated_base);
          causes[index] = snv.cause;

          ++index;
        }
      }
    }
  }
}

Rcpp::List PhylogeneticForest::get_sampled_cell_SNVs() const
{
  size_t num_of_mutations = count_SNVs(get_leaves_mutations());

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), chr_pos(num_of_mutations),
                alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), contexts(num_of_mutations),
                  mutated_bases(num_of_mutations), causes(num_of_mutations);

  size_t index{0};
  for (const auto& [cell_id, cell_mutations]: get_leaves_mutations()) {
    fill_lists(cell_mutations, cell_ids, chr_names, chr_pos, alleles,
               contexts, mutated_bases, causes, index);
  }

  return DataFrame::create(_["cell_id"]=cell_ids, _["chromosome"]=chr_names,
                           _["chr_pos"]=chr_pos, _["allele"]=alleles, 
                           _["context"]=contexts, _["mutated_base"]=mutated_bases,
                           _["cause"]=causes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_SNVs(const Races::Mutants::CellId& cell_id) const
{
  auto mutation_it = get_leaves_mutations().find(cell_id);

  if (mutation_it == get_leaves_mutations().end()) {
    throw std::domain_error("Cell \""+std::to_string(cell_id)+"\" is not a leaf");
  }

  const auto& cell_mutations = mutation_it->second;

  size_t num_of_mutations = count_SNVs(cell_mutations);

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), chr_pos(num_of_mutations),
                alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), contexts(num_of_mutations),
                  mutated_bases(num_of_mutations), causes(num_of_mutations);
  
  size_t index{0};
  fill_lists(cell_mutations, cell_ids, chr_names, chr_pos, alleles,
             contexts, mutated_bases, causes, index);

  return DataFrame::create(_["cell_id"]=cell_ids, _["chromosome"]=chr_names,
                           _["chr_pos"]=chr_pos, _["allele"]=alleles, 
                           _["context"]=contexts, _["mutated_base"]=mutated_bases,
                           _["cause"]=causes);
}

Races::Mutants::CellId PhylogeneticForest::get_first_occurrence(const SNV& snv) const
{
  auto first_cell_it = get_SNV_first_cell().find(snv);

  if (first_cell_it == get_SNV_first_cell().end()) {
    std::ostringstream oss;

    oss << "The mutation " << snv << " does not occurs in the sampled cells.";

    throw std::domain_error(oss.str());
  }

  return first_cell_it->second;
}

Races::Mutants::CellId PhylogeneticForest::get_first_occurrence(const CNA& cna) const
{
  auto first_cell_it = get_CNA_first_cell().find(cna);

  if (first_cell_it == get_CNA_first_cell().end()) {
    std::ostringstream oss;

    oss << "The mutation " << cna << " does not occurs in the sampled cells.";

    throw std::domain_error(oss.str());
  }

  return first_cell_it->second;
}

void PhylogeneticForest::save(const std::string& filename) const
{
  Races::Archive::Binary::Out out_archive(filename);

  Races::Mutations::PhylogeneticForest::save(out_archive);
}

PhylogeneticForest PhylogeneticForest::load(const std::string& filename)
{
  PhylogeneticForest forest;

  Races::Archive::Binary::In in_archive(filename);

  static_cast<Races::Mutations::PhylogeneticForest&>(forest) = Races::Mutations::PhylogeneticForest::load(in_archive);

  return forest;
}

void PhylogeneticForest::show() const
{
  using namespace Rcpp;

  size_t num_of_leaves{0};
  for (const auto& sample: get_samples()) {
    num_of_leaves += sample.get_cell_ids().size();
  }

  Rcout << "PhylogeneticForest" << std::endl 
        << "  # of trees: " << get_roots().size() << std::endl 
        << "  # of nodes: " << num_of_nodes() << std::endl 
        << "  # of leaves: " << num_of_leaves << std::endl 
        << "  samples: {";

  std::string sep = "";
  for (const auto& sample: get_samples()) {
    Rcout << sep << "\"" << sample.get_name() << "\"";
    sep = ", ";
  }

  Rcout << "}" << std::endl << std::endl
        << "  # of emerged SNVs: " << get_SNV_first_cell().size() << std::endl
        << "  # of emerged CNAs: " << get_CNA_first_cell().size() << std::endl
        << std::endl;
}

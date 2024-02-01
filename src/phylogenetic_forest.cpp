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

#include "phylogenetic_forest.hpp"

PhylogeneticForest::PhylogeneticForest():
  Races::Mutations::PhylogeneticForest()
{}

PhylogeneticForest::PhylogeneticForest(const Races::Mutations::PhylogeneticForest& orig,
                                       const std::filesystem::path& reference_path,
                                       const std::map<Races::Time, Races::Mutations::Exposure>& timed_exposures):
   Races::Mutations::PhylogeneticForest(orig), reference_path(reference_path), timed_exposures(timed_exposures)
{}

PhylogeneticForest::PhylogeneticForest(Races::Mutations::PhylogeneticForest&& orig,
                                       const std::filesystem::path& reference_path,
                                       const std::map<Races::Time, Races::Mutations::Exposure>& timed_exposures):
   Races::Mutations::PhylogeneticForest(std::move(orig)), reference_path(reference_path), timed_exposures(timed_exposures)
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

size_t count_CNAs(const Races::Mutations::CellGenomeMutations& cell_mutations)
{
  size_t counter{0};
  for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {
    counter += chromosome.get_CNAs().size();
  }

  return counter;
}

size_t count_CNAs(const std::map<Races::Mutants::CellId, Races::Mutations::CellGenomeMutations>& genome_mutations)
{
  size_t counter{0};
  for (const auto& [cell_id, cell_mutations]: genome_mutations) {
    counter += count_CNAs(cell_mutations);
  }

  return counter;
}

void fill_SNV_lists(const Races::Mutations::CellGenomeMutations& cell_mutations,
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

void fill_CNA_lists(const Races::Mutations::CellGenomeMutations& cell_mutations,
                    Rcpp::IntegerVector& cell_ids, Rcpp::CharacterVector& chr_names, 
                    Rcpp::IntegerVector& CNA_begins, Rcpp::IntegerVector& CNA_ends,
                    Rcpp::IntegerVector& src_alleles, Rcpp::IntegerVector& dst_alleles,
                    Rcpp::CharacterVector& types, size_t& index)
{
  using namespace Races::Mutations;

  for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {
    for (const auto& CNA: chromosome.get_CNAs()) {
      cell_ids[index] = cell_mutations.get_id();
      chr_names[index] = GenomicPosition::chrtos(chr_id);
      CNA_begins[index] = CNA.region.get_initial_position();
      CNA_ends[index] = CNA.region.get_final_position();
      bool is_amp = CNA.type == CopyNumberAlteration::Type::AMPLIFICATION;
      src_alleles[index] = (is_amp? CNA.source: NA_INTEGER);
      dst_alleles[index] = CNA.dest;

      types[index] = (is_amp?"A":"D");

      ++index;
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
    fill_SNV_lists(cell_mutations, cell_ids, chr_names, chr_pos, alleles,
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
  fill_SNV_lists(cell_mutations, cell_ids, chr_names, chr_pos, alleles,
                 contexts, mutated_bases, causes, index);

  return DataFrame::create(_["cell_id"]=cell_ids, _["chromosome"]=chr_names,
                           _["chr_pos"]=chr_pos, _["allele"]=alleles, 
                           _["context"]=contexts, _["mutated_base"]=mutated_bases,
                           _["cause"]=causes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_CNAs() const
{
  size_t num_of_mutations = count_CNAs(get_leaves_mutations());

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), CNA_begins(num_of_mutations),
                CNA_ends(num_of_mutations),
                src_alleles(num_of_mutations), dst_alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), types(num_of_mutations);

  size_t index{0};
  for (const auto& [cell_id, cell_mutations]: get_leaves_mutations()) {
    fill_CNA_lists(cell_mutations, cell_ids, chr_names, CNA_begins, CNA_ends,
                  src_alleles, dst_alleles, types, index);
  }

  return DataFrame::create(_["cell_id"]=cell_ids, _["type"]=types, 
                           _["chromosome"]=chr_names,
                           _["begin"]=CNA_begins, _["end"]=CNA_ends, 
                           _["allele"]=dst_alleles, _["src allele"]=src_alleles);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_CNAs(const Races::Mutants::CellId& cell_id) const
{
  auto mutation_it = get_leaves_mutations().find(cell_id);

  if (mutation_it == get_leaves_mutations().end()) {
    throw std::domain_error("Cell \""+std::to_string(cell_id)+"\" is not a leaf");
  }

  const auto& cell_mutations = mutation_it->second;
  size_t num_of_mutations = count_CNAs(cell_mutations);

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), CNA_begins(num_of_mutations),
                CNA_ends(num_of_mutations),
                src_alleles(num_of_mutations), dst_alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), types(num_of_mutations);
  
  size_t index{0};
  fill_CNA_lists(cell_mutations, cell_ids, chr_names, CNA_begins, CNA_ends,
                 src_alleles, dst_alleles, types, index);

  return DataFrame::create(_["cell_id"]=cell_ids, _["type"]=types, 
                           _["chromosome"]=chr_names,
                           _["begin"]=CNA_begins, _["end"]=CNA_ends, 
                           _["allele"]=dst_alleles, _["src allele"]=src_alleles);
}

template<typename MUTATION_TYPE, typename R_MUTATION> 
Rcpp::List get_first_occurrence(const std::map<MUTATION_TYPE, std::set<Races::Mutants::CellId>>& mutation_first_cells,
                                const R_MUTATION& mutation)
{
  auto first_cell_it = mutation_first_cells.find(mutation);

  if (first_cell_it == mutation_first_cells.end()) {
    std::ostringstream oss;

    oss << "The mutation " << mutation << " does not occurs in the sampled cells.";

    throw std::domain_error(oss.str());
  }

  const auto& cell_ids = first_cell_it->second;

  Rcpp::List R_cell_ids(cell_ids.size());

  auto ids_it = cell_ids.begin();
  for (size_t i=0; i<cell_ids.size(); ++i, ++ids_it) {
    R_cell_ids[i] = *ids_it;
  }
  return R_cell_ids;
}

Rcpp::List PhylogeneticForest::get_first_occurrence(const SEXP& mutation) const
{
  if ( TYPEOF(mutation) == S4SXP ) {

    Rcpp::S4 s4obj( mutation );
    if ( s4obj.is("Rcpp_SNV" ) ) {
      Rcpp::Environment env( s4obj );
      Rcpp::XPtr<SNV> snv_ptr( env.get(".pointer") );

      return ::get_first_occurrence(get_SNV_first_cells(), *snv_ptr);
    }

    if ( s4obj.is("Rcpp_CNA" ) ) {
      Rcpp::Environment env( s4obj );
      Rcpp::XPtr<CNA> cna_ptr( env.get(".pointer") );

    std::cout << *cna_ptr << std::endl;

      return ::get_first_occurrence(get_CNA_first_cells(), *cna_ptr);
    }
  }

  ::Rf_error("This methods exclusively accepts as the parameters "
             "SNV or CNA objects.");
}

Rcpp::List PhylogeneticForest::get_timed_exposures() const
{
  using namespace Rcpp;
  using namespace Races::Mutants;

  size_t dataframe_size{0};
  for (const auto& [time, exposure]: timed_exposures) {
    for (const auto& [SBS, prob]: exposure) {
      (void)SBS;

      ++dataframe_size;
    } 
  }

  NumericVector times(dataframe_size), probs(dataframe_size);
  CharacterVector SBSs(dataframe_size), types(dataframe_size);

  size_t index{0};
  for (const auto& [time, exposure]: timed_exposures) {
    for (const auto& [SBS, prob]: exposure) {
      times[index] = time;
      probs[index] = prob;
      SBSs[index] = SBS;
      types[index] = "SBS";
      ++index;
    }
  }

  return DataFrame::create(_["time"]=times, _["signature"]=SBSs,
                           _["exposure"]=probs, _["type"]=types);
}

void PhylogeneticForest::save(const std::string& filename) const
{
  Races::Archive::Binary::Out out_archive(filename);

  out_archive & reference_path
              & timed_exposures;

  out_archive.save(static_cast<const Races::Mutations::PhylogeneticForest&>(*this), "forest");
}

PhylogeneticForest PhylogeneticForest::load(const std::string& filename)
{
  Races::Archive::Binary::In in_archive(filename);

  std::string reference_path;

  in_archive & reference_path;

  PhylogeneticForest forest;
  forest.reference_path = reference_path;

  in_archive & forest.timed_exposures;

  in_archive.load(static_cast<Races::Mutations::PhylogeneticForest&>(forest), "forest");

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
        << "  # of emerged SNVs: " << get_SNV_first_cells().size() << std::endl
        << "  # of emerged CNAs: " << get_CNA_first_cells().size() << std::endl
        << std::endl;
}

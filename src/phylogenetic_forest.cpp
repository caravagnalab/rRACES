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

#include <utils.hpp>

#include "phylogenetic_forest.hpp"
#include "sid.hpp"
#include "cna.hpp"
#include "simulation.hpp"

#include "utility.hpp"

PhylogeneticForest::PhylogeneticForest():
  RACES::Mutations::PhylogeneticForest()
{}

PhylogeneticForest::PhylogeneticForest(const RACES::Mutations::PhylogeneticForest& orig,
                                       const std::string& germline_subject,
                                       const std::filesystem::path reference_path,
                                       const TimedMutationalExposure& timed_SBS_exposures,
                                       const TimedMutationalExposure& timed_indel_exposures):
    RACES::Mutations::PhylogeneticForest(orig), germline_subject(germline_subject),
    reference_path(reference_path)
{
    using namespace RACES::Mutations;
    timed_exposures[MutationType::Type::SBS] = timed_SBS_exposures;
    timed_exposures[MutationType::Type::INDEL] = timed_indel_exposures;
}

PhylogeneticForest::PhylogeneticForest(RACES::Mutations::PhylogeneticForest&& orig,
                                       const std::string& germline_subject,
                                       const std::filesystem::path reference_path,
                                       const TimedMutationalExposure& timed_SBS_exposures,
                                       const TimedMutationalExposure& timed_indel_exposures):
   RACES::Mutations::PhylogeneticForest(std::move(orig)), germline_subject(germline_subject),
   reference_path(reference_path)
{
    using namespace RACES::Mutations;
    timed_exposures[MutationType::Type::SBS] = timed_SBS_exposures;
    timed_exposures[MutationType::Type::INDEL] = timed_indel_exposures;
}

Rcpp::List PhylogeneticForest::get_samples_info() const
{
return SpatialSimulation::get_samples_info(get_samples());
}

PhylogeneticForest PhylogeneticForest::get_subforest_for(const std::vector<std::string>& sample_names) const
{
  PhylogeneticForest forest;

  forest.reference_path = reference_path;
  forest.timed_exposures = timed_exposures;

  static_cast<RACES::Mutations::PhylogeneticForest&>(forest) = RACES::Mutations::PhylogeneticForest::get_subforest_for(sample_names);

  return forest;
}

size_t count_mutations(const RACES::Mutations::GenomeMutations& mutations)
{
  size_t counter{0};
  for (const auto& [chr_id, chromosome]: mutations.get_chromosomes()) {
    for (const auto& [allele_id, allele]: chromosome.get_alleles()) {
      for (const auto& [fragment_pos, fragment]: allele.get_fragments()) {
          counter += fragment.get_mutations().size();
      }
    }
  }

  return counter;
}

size_t count_mutations(const std::map<RACES::Mutants::CellId, std::shared_ptr<RACES::Mutations::CellGenomeMutations>>& genome_mutations)
{
  size_t counter{0};
  for (const auto& [cell_id, mutations_ptr]: genome_mutations) {
    counter += count_mutations(*mutations_ptr);
  }

  return counter;
}

size_t count_CNAs(const RACES::Mutations::GenomeMutations& mutations)
{
  size_t counter{0};
  for (const auto& [chr_id, chromosome]: mutations.get_chromosomes()) {
    counter += chromosome.get_CNAs().size();
  }

  return counter;
}

size_t count_CNAs(const std::map<RACES::Mutants::CellId, std::shared_ptr<RACES::Mutations::CellGenomeMutations>>& genome_mutations)
{
  size_t counter{0};
  for (const auto& [cell_id, mutations_ptr]: genome_mutations) {
    counter += count_CNAs(*mutations_ptr);
  }

  return counter;
}

void fill_mutation_lists(const RACES::Mutations::CellGenomeMutations& cell_mutations,
                         Rcpp::IntegerVector& cell_ids, Rcpp::CharacterVector& chr_names,
                         Rcpp::IntegerVector& chr_pos, Rcpp::IntegerVector& alleles,
                         Rcpp::CharacterVector& ref, Rcpp::CharacterVector& alt,
                         Rcpp::CharacterVector& types, Rcpp::CharacterVector& causes,
                         Rcpp::CharacterVector& classes, size_t& index)
{
  using namespace RACES::Mutations;

  for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {
    for (const auto& [allele_id, allele]: chromosome.get_alleles()) {
      for (const auto& [fragment_pos, fragment]: allele.get_fragments()) {
        for (const auto& [mutation_pos, mutation]: fragment.get_mutations()) {
          cell_ids[index] = cell_mutations.get_id();
          chr_names[index] = GenomicPosition::chrtos(chr_id);
          chr_pos[index] = mutation.position;
          alleles[index] = allele_id;
          ref[index] = mutation.ref;
          alt[index] = mutation.alt;
          types[index] = (mutation.is_SBS()?"SNV":"indel");
          causes[index] = mutation.cause;
          classes[index] = mutation.get_nature_description();

          ++index;
        }
      }
    }
  }
}

void fill_CNA_lists(const RACES::Mutations::CellGenomeMutations& cell_mutations,
                    Rcpp::IntegerVector& cell_ids, Rcpp::CharacterVector& chr_names,
                    Rcpp::IntegerVector& CNA_begins, Rcpp::IntegerVector& CNA_ends,
                    Rcpp::IntegerVector& src_alleles, Rcpp::IntegerVector& dst_alleles,
                    Rcpp::CharacterVector& types, Rcpp::CharacterVector& classes,
                    size_t& index)
{
  for (const auto& [chr_id, chromosome]: cell_mutations.get_chromosomes()) {
    for (const auto& cna: chromosome.get_CNAs()) {
      cell_ids[index] = cell_mutations.get_id();
      chr_names[index] = RACES::Mutations::GenomicPosition::chrtos(chr_id);
      CNA_begins[index] = cna.get_initial_position();
      CNA_ends[index] = cna.get_final_position();
      bool is_amp = cna.type == RACES::Mutations::CNA::Type::AMPLIFICATION;
      src_alleles[index] = (is_amp? cna.source: NA_INTEGER);
      dst_alleles[index] = cna.dest;
      classes[index] = cna.get_nature_description();

      types[index] = (is_amp?"A":"D");

      ++index;
    }
  }
}

Rcpp::List PhylogeneticForest::get_absolute_chromosome_positions() const
{
  const auto& germline = get_germline_mutations();

  auto abspos = germline.get_absolute_chromosome_positions();

  using namespace Rcpp;

  NumericVector lengths(abspos.size()), froms(abspos.size()),
                tos(abspos.size());
  CharacterVector names(abspos.size());

  size_t index{0};
  size_t pos{1};
  for (const auto& [chr_id, from]: abspos) {
    names[index] = RACES::Mutations::GenomicPosition::chrtos(chr_id);
    auto chr_length = germline.get_chromosome(chr_id).size();
    lengths[index] = chr_length;
    froms[index] = pos;
    pos += chr_length;
    tos[index] = (pos - 1);
    ++index;
  }

  return DataFrame::create(_["chr"]=names,
                           _["length"]=lengths, _["from"]=froms,
                           _["to"]=tos);
}

Rcpp::List PhylogeneticForest::get_germline_SIDs() const
{
  const auto& germline = get_germline_mutations();
  size_t num_of_mutations{0};
  for (const auto& [chr_id, chr] : germline.get_chromosomes()) {
    for (const auto& [allele_id, allele] : chr.get_alleles()) {
      for (const auto& [f_pos, fragment] : allele.get_fragments()) {
        num_of_mutations += fragment.get_mutations().size();
      }
    }
  }

  using namespace Rcpp;

  IntegerVector chr_pos(num_of_mutations), alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
                  alt(num_of_mutations), types(num_of_mutations),
                  causes(num_of_mutations), classes(num_of_mutations);

  size_t index{0};
  for (const auto& [chr_id, chr] : germline.get_chromosomes()) {
    for (const auto& [allele_id, allele] : chr.get_alleles()) {
      for (const auto& [f_pos, fragment] : allele.get_fragments()) {
        for (const auto& [mutation_pos, mutation] : fragment.get_mutations()) {
          chr_names[index] = RACES::Mutations::GenomicPosition::chrtos(chr_id);
          chr_pos[index] = mutation.position;
          alleles[index] = allele_id;
          ref[index] = mutation.ref;
          alt[index] = mutation.alt;
          types[index] = (mutation.is_SBS()?"SNV":"indel");
          causes[index] = mutation.cause;
          classes[index] = mutation.get_nature_description();

          ++index;
        }
      }
    }
  }

  return DataFrame::create(_["chr"]=chr_names,
                           _["chr_pos"]=chr_pos, _["allele"]=alleles,
                           _["ref"]=ref, _["alt"]=alt, _["type"]=types,
                           _["cause"]=causes, _["class"]=classes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_SIDs() const
{
  size_t num_of_mutations = count_mutations(get_leaves_mutations());

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), chr_pos(num_of_mutations),
                alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
                  alt(num_of_mutations), types(num_of_mutations),
                  causes(num_of_mutations), classes(num_of_mutations);

  size_t index{0};
  for (const auto& [cell_id, mutations_ptr]: get_leaves_mutations()) {
    fill_mutation_lists(*mutations_ptr, cell_ids, chr_names, chr_pos,
                        alleles, ref, alt, types, causes, classes, index);
  }

  return DataFrame::create(_["cell_id"]=cell_ids, _["chr"]=chr_names,
                           _["chr_pos"]=chr_pos, _["allele"]=alleles,
                           _["ref"]=ref, _["alt"]=alt, _["type"]=types,
                           _["cause"]=causes, _["class"]=classes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_SIDs(const RACES::Mutants::CellId& cell_id) const
{
  auto mutation_it = get_leaves_mutations().find(cell_id);

  if (mutation_it == get_leaves_mutations().end()) {
    throw std::domain_error("Cell \""+std::to_string(cell_id)+"\" is not a leaf");
  }

  return get_SID_dataframe(*(mutation_it->second));
}

Rcpp::List PhylogeneticForest::get_SID_dataframe(const RACES::Mutations::CellGenomeMutations& cell_mutations)
{
  size_t num_of_mutations = count_mutations(cell_mutations);

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), chr_pos(num_of_mutations),
                alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
                  alt(num_of_mutations), types(num_of_mutations),
                  causes(num_of_mutations), classes(num_of_mutations);

  size_t index{0};
  fill_mutation_lists(cell_mutations, cell_ids, chr_names, chr_pos, alleles,
                      ref, alt, types, causes, classes, index);

  return DataFrame::create(_["cell_id"]=cell_ids, _["chr"]=chr_names,
                           _["chr_pos"]=chr_pos, _["allele"]=alleles,
                           _["ref"]=ref, _["alt"]=alt, _["type"]=types,
                           _["cause"]=causes, _["class"]=classes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_CNAs() const
{
  size_t num_of_mutations = count_CNAs(get_leaves_mutations());

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), CNA_begins(num_of_mutations),
                CNA_ends(num_of_mutations),
                src_alleles(num_of_mutations), dst_alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), types(num_of_mutations),
                  classes(num_of_mutations);

  size_t index{0};
  for (const auto& [cell_id, mutations_ptr]: get_leaves_mutations()) {
    fill_CNA_lists(*mutations_ptr, cell_ids, chr_names, CNA_begins, CNA_ends,
                   src_alleles, dst_alleles, types, classes, index);
  }

  return DataFrame::create(_["cell_id"]=cell_ids, _["type"]=types,
                           _["chr"]=chr_names,
                           _["begin"]=CNA_begins, _["end"]=CNA_ends,
                           _["allele"]=dst_alleles, _["src allele"]=src_alleles,
                           _["class"]=classes);
}

Rcpp::List PhylogeneticForest::get_sampled_cell_CNAs(const RACES::Mutants::CellId& cell_id) const
{
  auto mutation_it = get_leaves_mutations().find(cell_id);

  if (mutation_it == get_leaves_mutations().end()) {
    throw std::domain_error("Cell \""+std::to_string(cell_id)+"\" is not a leaf");
  }

  const auto& cell_mutations = *(mutation_it->second);
  size_t num_of_mutations = count_CNAs(cell_mutations);

  using namespace Rcpp;

  IntegerVector cell_ids(num_of_mutations), CNA_begins(num_of_mutations),
                CNA_ends(num_of_mutations),
                src_alleles(num_of_mutations), dst_alleles(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), types(num_of_mutations),
                  classes(num_of_mutations);

  size_t index{0};
  fill_CNA_lists(cell_mutations, cell_ids, chr_names, CNA_begins, CNA_ends,
                 src_alleles, dst_alleles, types, classes, index);

  return DataFrame::create(_["cell_id"]=cell_ids, _["type"]=types,
                           _["chr"]=chr_names,
                           _["begin"]=CNA_begins, _["end"]=CNA_ends,
                           _["allele"]=dst_alleles, _["src allele"]=src_alleles,
                           _["class"]=classes);
}

template<typename MUTATION_TYPE, typename R_MUTATION>
Rcpp::List get_first_occurrence(const std::map<MUTATION_TYPE, std::set<RACES::Mutants::CellId>>& mutation_first_cells,
                                const RACES::Mutations::GenomeMutations& germline,
                                const R_MUTATION& mutation)
{
  auto first_cell_it = mutation_first_cells.find(mutation);

  if (first_cell_it == mutation_first_cells.end()) {
    std::ostringstream oss;

    if constexpr (std::is_base_of_v<MUTATION_TYPE, RACES::Mutations::SID>) {
      if (germline.includes(mutation)) {
        Rcpp::List R_cell_ids(1);

        R_cell_ids[0] = NA_INTEGER;
        return R_cell_ids;
      } else {
        oss << "The mutation " << mutation << " does not occurs in the sampled cells.";
      }
    } else {
      oss << "The mutation " << mutation << " does not occurs in the sampled cells.";
    }

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
    if ( s4obj.is("Rcpp_Mutation" ) ) {
      Rcpp::Environment env( s4obj );
      Rcpp::XPtr<SID> snv_ptr( env.get(".pointer") );

      return ::get_first_occurrence(get_mutation_first_cells(), get_germline_mutations(),
                                    *snv_ptr);
    }

    if ( s4obj.is("Rcpp_CNA" ) ) {
      Rcpp::Environment env( s4obj );
      Rcpp::XPtr<CNA> cna_ptr( env.get(".pointer") );

      return ::get_first_occurrence(get_CNA_first_cells(), get_germline_mutations(),
                                    *cna_ptr);
    }
  }

  ::Rf_error("This methods exclusively accepts as the parameters "
             "SNV, Indel, or CNA objects.");
}

void fill_timed_exposures(const std::map<RACES::Time, RACES::Mutations::MutationalExposure>& mutation_timed_exposures,
                          const std::string& mutation_type_name,
                          size_t& index, Rcpp::NumericVector& times, Rcpp::NumericVector& probs,
                          Rcpp::CharacterVector& sig_names, Rcpp::CharacterVector& types)
{
  for (const auto& [time, exposure]: mutation_timed_exposures) {
    for (const auto& [sign_name, prob]: exposure) {
      times[index] = time;
      probs[index] = prob;
      sig_names[index] = sign_name;
      types[index] = mutation_type_name;
      ++index;
    }
  }
}


Rcpp::List PhylogeneticForest::get_timed_exposures() const
{
  using namespace Rcpp;
  using namespace RACES::Mutants;
  using namespace RACES::Mutations;

  size_t dataframe_size{0};
  for (const auto& [type, mutation_timed_exposures] : timed_exposures) {
    (void)type;
    for (const auto& [time, exposure]: mutation_timed_exposures) {
      for (const auto& [sign_name, prob]: exposure) {
        (void)sign_name;
        ++dataframe_size;
      }
    }
  }

  NumericVector times(dataframe_size), probs(dataframe_size);
  CharacterVector sig_names(dataframe_size), types(dataframe_size);

  size_t index{0};
  fill_timed_exposures(timed_exposures.at(MutationType::Type::SBS), "SNV",
                       index, times, probs, sig_names, types);
  fill_timed_exposures(timed_exposures.at(MutationType::Type::INDEL), "indel",
                       index, times, probs, sig_names, types);

  return DataFrame::create(_["time"]=times, _["signature"]=sig_names,
                           _["exposure"]=probs, _["type"]=types);
}

void PhylogeneticForest::set_reference_path(const std::string reference_path)
{
  if (!std::filesystem::exists(reference_path)) {
    throw std::runtime_error("The reference genome file \""+ reference_path
                             + "\" does not exists.");
  }

  this->reference_path = std::filesystem::absolute(std::filesystem::path(reference_path));
}

void PhylogeneticForest::save(const std::string& filename) const
{
  RACES::Archive::Binary::Out out_archive(filename);

  auto ref_str = to_string(reference_path);

  out_archive & ref_str
              & timed_exposures;

  out_archive.save(static_cast<const RACES::Mutations::PhylogeneticForest&>(*this),
                   "forest", Rcpp::Rcout);
}

PhylogeneticForest PhylogeneticForest::load(const std::string& filename)
{
  RACES::Archive::Binary::In in_archive(filename);

  std::string reference_path;

  in_archive & reference_path;

  PhylogeneticForest forest;
  forest.reference_path = reference_path;

  in_archive & forest.timed_exposures;

  try {
    in_archive.load(static_cast<RACES::Mutations::PhylogeneticForest&>(forest),
                    "forest", Rcpp::Rcout);
  } catch (RACES::Archive::WrongFileFormatDescr& ex) {
    raise_error(ex, "phylogenetic forest");
  } catch (RACES::Archive::WrongFileFormatVersion& ex) {
    raise_error(ex, "phylogenetic forest");
  }

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
        << "  # of emerged SNVs and indels: " << get_mutation_first_cells().size() << std::endl
        << "  # of emerged CNAs: " << get_CNA_first_cells().size() << std::endl
        << std::endl;
}

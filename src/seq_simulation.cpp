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

#include "seq_simulation.hpp"
#include "sequencers.hpp"
#include "sampled_cell.hpp"

#include "utility.hpp"

std::string join(const std::set<std::string>& S, const char& sep=';')
{
  std::ostringstream oss;

  if (S.size()>0) {
    auto S_it = S.begin();
    oss << *S_it;
    while (++S_it != S.end()) {
      oss << sep << *S_it;
    }
  }

  return oss.str();
}

std::set<std::string> get_descriptions(const std::set<RACES::Mutations::Mutation::Nature>& nature_set)
{
  std::set<std::string> nature_strings;

  for (const auto& nature: nature_set) {
    nature_strings.insert(RACES::Mutations::Mutation::get_nature_description(nature));
  }

  return nature_strings;
}

void add_SNV_data(Rcpp::DataFrame& df,
                  const RACES::Mutations::SequencingSimulations::SampleStatistics& sample_statistics)
{
  using namespace Rcpp;
  using namespace RACES::Mutations;

  size_t num_of_mutations = sample_statistics.get_data().size();

  IntegerVector chr_pos(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), ref(num_of_mutations),
                  alt(num_of_mutations), causes(num_of_mutations),
                  classes(num_of_mutations);

  size_t index{0};
  for (const auto& [mutation, data] : sample_statistics.get_data()) {
    chr_names[index] = GenomicPosition::chrtos(mutation.chr_id);
    chr_pos[index] = mutation.position;
    ref[index] = mutation.ref;
    alt[index] = mutation.alt;

    auto full_causes = join(data.causes, ';');

    if (full_causes == "") {
      causes[index] = NA_STRING;
    } else {
      causes[index] = full_causes;
    }

    auto descr = get_descriptions(data.nature_set);
    classes[index] = join(descr, ';');

    ++index;
  }

  df.push_back(chr_names, "chr");
  df.push_back(chr_pos, "chr_pos");
  df.push_back(ref, "ref");
  df.push_back(alt, "alt");
  df.push_back(causes, "causes");
  df.push_back(classes, "classes");
}

void add_sample_statistics(Rcpp::DataFrame& df,
                           const RACES::Mutations::SequencingSimulations::SampleStatistics& sample_statistics)
{
  if (df.nrows()==0) {
    add_SNV_data(df, sample_statistics);
  }

  size_t num_of_mutations = sample_statistics.get_data().size();

  if (num_of_mutations != static_cast<size_t>(df.nrows())) {
    throw std::runtime_error("SeqSimResults are not canonical!!!");
  }

  using namespace Rcpp;
  using namespace RACES::Mutations;

  DoubleVector VAF(num_of_mutations);
  IntegerVector occurrences(num_of_mutations), coverages(num_of_mutations);

  size_t index{0};
  auto coverage_it = sample_statistics.get_coverage().begin();
  std::less<GenomicPosition> come_before;
  for (const auto& [mutation, mutation_data] : sample_statistics.get_data()) {
    occurrences[index] = mutation_data.num_of_occurrences;

    if (come_before(coverage_it->first, mutation)) {
      ++coverage_it;
    }

    coverages[index] = coverage_it->second;
    VAF[index] = static_cast<double>(mutation_data.num_of_occurrences)/coverage_it->second;

    ++index;
  }

  const auto& sample_name = sample_statistics.get_sample_name();

  df.push_back(occurrences, sample_name+".occurrences");
  df.push_back(coverages, sample_name+".coverage");
  df.push_back(VAF, sample_name+".VAF");
}

Rcpp::List get_result_dataframe(const RACES::Mutations::SequencingSimulations::SampleSetStatistics& sample_set_statistics)
{
  auto df = Rcpp::DataFrame::create();

  for (const auto& [sample_name, sample_stats] : sample_set_statistics) {
    add_sample_statistics(df, sample_stats);
  }

  return df;
}

void
split_by_labels(std::list<RACES::Mutations::SampleGenomeMutations>& FACS_samples,
                const Rcpp::Function& labelling_function,
                const RACES::Mutations::SampleGenomeMutations& sample_mutations,
                const PhylogeneticForest& forest)
{
    using namespace RACES::Mutants;
    using namespace RACES::Mutants::Evolutions;
    using namespace RACES::Mutations;

    std::map<std::string, SampleGenomeMutations*> labelled_samples;

    for (const auto& cell_mutations : sample_mutations.mutations) {
        auto node = Rcpp::wrap(SampledCell(forest, cell_mutations->get_id()));

        std::string label = Rcpp::as<std::string>(labelling_function(node));

        auto found = labelled_samples.find(label);

        if (found == labelled_samples.end()) {
            std::string new_name = sample_mutations.name;
            if (label != "") {
                new_name = new_name + "_" + label;
            }

            FACS_samples.emplace_back(new_name, sample_mutations.germline_mutations);

            FACS_samples.back().mutations.push_back(cell_mutations);

            labelled_samples.insert({label, &(FACS_samples.back())});
        } else {
            (found->second)->mutations.push_back(cell_mutations);
        }
    }
}

void
apply_FACS_labels(std::list<RACES::Mutations::SampleGenomeMutations>& sample_mutations_list,
                  const SEXP& labelling_function,
                  const PhylogeneticForest& forest)
{
    switch (TYPEOF(labelling_function)) {
        case NILSXP:
            break;
        case CLOSXP:
        {
            std::list<RACES::Mutations::SampleGenomeMutations> FACS_samples;
            Rcpp::Function l_function = Rcpp::as<Rcpp::Function>(labelling_function);

            for (const auto& sample_mutations : sample_mutations_list) {
                split_by_labels(FACS_samples, l_function, sample_mutations, forest);
            }

            std::swap(sample_mutations_list, FACS_samples);

            break;
        }
        default:
            throw std::domain_error("The FACs_labelling_function must be a function.");
    } 
}

std::set<RACES::Mutations::ChromosomeId>
get_genome_chromosome_ids(const std::list<RACES::Mutations::SampleGenomeMutations>& mutations_list)
{
    for (const auto& sample_mutations : mutations_list) {
        for (const auto& mutations : sample_mutations.mutations) {
            std::set<RACES::Mutations::ChromosomeId> ids;
            for (const auto& [chr_id, chr_mutations] : mutations->get_chromosomes()) {
                ids.insert(chr_id);
            }

            return ids;
        }
    }

    return {};
}

std::set<RACES::Mutations::ChromosomeId>
get_relevant_chr_set(std::list<RACES::Mutations::SampleGenomeMutations> mutations_list,
                     SEXP& chromosome_ids)
{
  using namespace Rcpp;
  using namespace RACES::Mutations;

  switch (TYPEOF(chromosome_ids)) {
    case NILSXP:
    {
        return get_genome_chromosome_ids(mutations_list);
    }
    case STRSXP:
    {
        std::set<ChromosomeId> chr_ids;

        CharacterVector chr_names{chromosome_ids};

        for (const auto& chr_name : chr_names) {
            chr_ids.insert(GenomicPosition::stochr(as<std::string>(chr_name)));
        }
        return chr_ids;
    }
    case VECSXP:
    {
        std::set<ChromosomeId> chr_ids;
        List chr_names = Rcpp::as<List>(chromosome_ids);

        size_t i{0};
        for (const auto& chr_name : chr_names) {
            ++i;
            if (TYPEOF(chr_name) != STRSXP) {
                throw std::domain_error("Expected a list of string: the "
                                        + ordtostr(i)
                                        + " element of the list is not "
                                        + "a string.");
            }

            Rcpp::CharacterVector name{chr_name};
            if (name.length()>1) {
                throw std::domain_error("Expected a list of string: the "
                                        + ordtostr(i)
                                        + " element of the list is not "
                                        + "a string.");
            }

            chr_ids.insert(GenomicPosition::stochr(as<std::string>(name)));
        }
        return chr_ids;
    }
    default:
        throw std::domain_error("Unsupported sequencer type");
  }
}

RACES::Mutations::SequencingSimulations::SampleSetStatistics
simulate_seq(RACES::Mutations::SequencingSimulations::ReadSimulator<>& simulator,
             SEXP& sequencer,
             std::list<RACES::Mutations::SampleGenomeMutations> mutations_list,
             const std::set<RACES::Mutations::ChromosomeId>& chromosome_ids,
             const double& coverage, const double purity,
             const std::string& base_name, std::ostream& progress_bar_stream)
{
  switch (TYPEOF(sequencer)) {
    case S4SXP:
    {
      Rcpp::S4 s4obj( sequencer );
      if ( s4obj.is("Rcpp_BasicIlluminaSequencer")) {
        Rcpp::Environment env( s4obj );

        Rcpp::XPtr<BasicIlluminaSequencer> sequencer_ptr( env.get(".pointer") );

        return simulator(*sequencer_ptr, mutations_list, chromosome_ids,
                         coverage, purity, base_name, progress_bar_stream);
      }
      if ( s4obj.is("Rcpp_ErrorlessIlluminaSequencer")) {
        Rcpp::Environment env( s4obj );

        Rcpp::XPtr<ErrorlessIlluminaSequencer> sequencer_ptr( env.get(".pointer") );

        return simulator(*sequencer_ptr, mutations_list, chromosome_ids,
                         coverage, purity, base_name, progress_bar_stream);
      }
    }
    case NILSXP:
    {
        ErrorlessIlluminaSequencer sequencer;

        return simulator(sequencer, mutations_list, chromosome_ids, coverage,
                         purity, base_name, progress_bar_stream);
    }
    default:
        throw std::domain_error("Unsupported sequencer type");
  }
}

Rcpp::List simulate_seq(const PhylogeneticForest& forest, SEXP& sequencer,
                        SEXP& chromosome_ids, const double& coverage,
                        const int& read_size, const int& insert_size,
                        const std::string& output_dir, const bool& write_SAM,
                        const bool& update_SAM_dir,
                        const SEXP& FACS_labelling_function,
                        const double& purity, const bool& with_normal_sample,
                        const int& rnd_seed)
{
  using namespace RACES::Mutations::SequencingSimulations;

  ReadSimulator<> simulator;

  if (!std::filesystem::exists(forest.get_reference_path())) {
    throw std::runtime_error("The reference genome file \""
                             + to_string(forest.get_reference_path())
                             + "\" does not exists anymore. Please, re-build "
                             + "the mutation engine.");
  }

  std::filesystem::path output_path = output_dir;

  bool remove_output_path = false;

  if (!write_SAM) {
    remove_output_path = true;
    output_path = get_tmp_dir_path(output_dir);
  }

  ReadSimulator<>::Mode SAM_mode = ReadSimulator<>::Mode::CREATE;

  if (update_SAM_dir) {
    SAM_mode = ReadSimulator<>::Mode::UPDATE;
  }

  if (insert_size==0) {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                SAM_mode, rnd_seed);
  } else {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                insert_size, SAM_mode, rnd_seed);
  }

  simulator.enable_SAM_writing(write_SAM);

  auto mutations_list = forest.get_sample_mutations_list();

  apply_FACS_labels(mutations_list, FACS_labelling_function, forest);

  if (with_normal_sample) {
    const auto& germline = forest.get_germline_mutations();
    mutations_list.emplace_back("normal_sample", germline);
    auto germline_structure_ptr = std::make_shared<RACES::Mutations::CellGenomeMutations>(germline.duplicate_structure());
    mutations_list.back().mutations.push_back(germline_structure_ptr);
  }

  const auto chr_ids = get_relevant_chr_set(mutations_list, chromosome_ids);

  auto result = simulate_seq(simulator, sequencer, mutations_list, chr_ids, coverage,
                             purity, "chr_", Rcpp::Rcout);

  if (remove_output_path) {
    std::filesystem::remove_all(output_path);
  }

  return get_result_dataframe(result);
}

Rcpp::List simulate_normal_seq(const PhylogeneticForest& forest, SEXP& sequencer,
                               SEXP& chromosome_ids, const double& coverage,
                               const int& read_size, const int& insert_size,
                               const std::string& output_dir, const bool& write_SAM,
                               const bool& update_SAM_dir, const int& rnd_seed)
{
  using namespace RACES::Mutations::SequencingSimulations;

  ReadSimulator<> simulator;

  if (!std::filesystem::exists(forest.get_reference_path())) {
    throw std::runtime_error("The reference genome file \""
                             + to_string(forest.get_reference_path())
                             + "\" does not exists anymore. Please, re-build "
                             + "the mutation engine.");
  }

  std::filesystem::path output_path = output_dir;

  bool remove_output_path = false;

  if (!write_SAM) {
    remove_output_path = true;
    output_path = get_tmp_dir_path(output_dir);
  }

  ReadSimulator<>::Mode SAM_mode = ReadSimulator<>::Mode::CREATE;

  if (update_SAM_dir) {
    SAM_mode = ReadSimulator<>::Mode::UPDATE;
  }

  if (insert_size==0) {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                SAM_mode, rnd_seed);
  } else {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                insert_size, SAM_mode, rnd_seed);
  }

  simulator.enable_SAM_writing(write_SAM);

  const auto& germline = forest.get_germline_mutations();

  std::list<RACES::Mutations::SampleGenomeMutations> mutations_list;
  mutations_list.emplace_back("normal_sample", germline);
  auto germline_structure_ptr = std::make_shared<RACES::Mutations::CellGenomeMutations>(germline.duplicate_structure());
  mutations_list.front().mutations.push_back(germline_structure_ptr);

  const auto chr_ids = get_relevant_chr_set(mutations_list, chromosome_ids);

  auto result = simulate_seq(simulator, sequencer, mutations_list,  chr_ids, coverage,
                             0, "chr_", Rcpp::Rcout);

  if (remove_output_path) {
    std::filesystem::remove_all(output_path);
  }

  return get_result_dataframe(result);
}
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

std::set<std::string> get_descriptions(const std::set<Races::Mutations::Mutation::Nature>& nature_set)
{
  std::set<std::string> nature_strings;

  for (const auto& nature: nature_set) {
    nature_strings.insert(Races::Mutations::Mutation::get_nature_description(nature));
  }

  return nature_strings;
}

void add_SNV_data(Rcpp::DataFrame& df,
                  const Races::Mutations::SequencingSimulations::SampleStatistics& sample_statistics)
{
  using namespace Rcpp;
  using namespace Races::Mutations;

  size_t num_of_mutations = sample_statistics.get_SNV_data().size();

  IntegerVector chr_pos(num_of_mutations);
  CharacterVector chr_names(num_of_mutations), ref_bases(num_of_mutations),
                  alt_bases(num_of_mutations), causes(num_of_mutations),
                  classes(num_of_mutations);

  size_t index{0};
  for (const auto& [snv, data] : sample_statistics.get_SNV_data()) {
    chr_names[index] = GenomicPosition::chrtos(snv.chr_id);
    chr_pos[index] = snv.position;
    ref_bases[index] = std::string(1,snv.ref_base);
    alt_bases[index] = std::string(1,snv.alt_base);
  
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

  df.push_back(chr_names, "chromosome");
  df.push_back(chr_pos, "chr_pos");
  df.push_back(ref_bases, "ref");
  df.push_back(alt_bases, "alt");
  df.push_back(causes, "causes");
  df.push_back(classes, "classes");
}

void add_sample_statistics(Rcpp::DataFrame& df,
                           const Races::Mutations::SequencingSimulations::SampleStatistics& sample_statistics)
{
  if (df.nrows()==0) {
    add_SNV_data(df, sample_statistics);
  }

  size_t num_of_mutations = sample_statistics.get_SNV_data().size();

  if (num_of_mutations != static_cast<size_t>(df.nrows())) {
    throw std::runtime_error("SeqSimResults are not canonical!!!");
  }

  using namespace Rcpp;
  using namespace Races::Mutations;

  DoubleVector VAF(num_of_mutations);
  IntegerVector occurrences(num_of_mutations), coverages(num_of_mutations);

  size_t index{0};
  auto coverage_it = sample_statistics.get_SNV_coverage().begin();
  std::less<GenomicPosition> come_before; 
  for (const auto& [snv, snv_data] : sample_statistics.get_SNV_data()) {
    occurrences[index] = snv_data.num_of_occurrences;
  
    if (come_before(coverage_it->first, snv)) {
      ++coverage_it;
    }
  
    coverages[index] = coverage_it->second;
    VAF[index] = static_cast<double>(snv_data.num_of_occurrences)/coverage_it->second;
  
    ++index;
  }

  const auto& sample_name = sample_statistics.get_sample_name();

  df.push_back(occurrences, sample_name+".occurrences");
  df.push_back(coverages, sample_name+".coverage");
  df.push_back(VAF, sample_name+".VAF");
}

Rcpp::List get_result_dataframe(const Races::Mutations::SequencingSimulations::SampleSetStatistics& sample_set_statistics)
{
  auto df = Rcpp::DataFrame::create();

  for (const auto& [sample_name, sample_stats] : sample_set_statistics) {
    add_sample_statistics(df, sample_stats);
  }

  return df;
}

void
split_by_epigenetic_status(std::list<Races::Mutations::SampleGenomeMutations>& FACS_samples,
                           const Races::Mutations::SampleGenomeMutations& sample_mutations,
                           std::map<Races::Mutants::SpeciesId, std::string> methylation_map)
{
    using namespace Races::Mutants;
    using namespace Races::Mutants::Evolutions;
    using namespace Races::Mutations;

    std::map<SpeciesId, SampleGenomeMutations*> meth_samples;

    for (const auto& cell_mutations : sample_mutations.mutations) {
        auto found = meth_samples.find(cell_mutations->get_species_id());

        if (found == meth_samples.end()) {
            auto new_name = sample_mutations.name+"_"+
                                methylation_map.at(cell_mutations->get_species_id());

            FACS_samples.emplace_back(new_name, sample_mutations.germline_mutations);

            FACS_samples.back().mutations.push_back(cell_mutations);

            meth_samples.insert({cell_mutations->get_species_id(), &(FACS_samples.back())});
        } else {
            (found->second)->mutations.push_back(cell_mutations);
        }
    }
}

std::list<Races::Mutations::SampleGenomeMutations>
split_by_epigenetic_status(const std::list<Races::Mutations::SampleGenomeMutations>& sample_mutations_list,
                           const PhylogeneticForest& forest)
{
    using namespace Races::Mutants;
    using namespace Races::Mutants::Evolutions;

    std::map<SpeciesId, std::string> methylation_map;

    for (const auto& [species_id, species_data] : forest.get_species_data()) {
      if (MutantProperties::signature_to_string(species_data.signature) == "+") {
        methylation_map[species_id] = "P";
      } else {
        methylation_map[species_id] = "N";
      }
    }

    std::list<Races::Mutations::SampleGenomeMutations> FACS_samples;

    for (const auto& sample_mutations : sample_mutations_list) {
        split_by_epigenetic_status(FACS_samples, sample_mutations, methylation_map);
    }

    return FACS_samples;
}

Rcpp::List simulate_seq(const PhylogeneticForest& forest, const double& coverage, 
                        const int& read_size, const int& insert_size,
                        const std::string& output_dir, const bool& write_SAM,
                        const bool& FACS, const double& purity, const int& rnd_seed)
{
  using namespace Races::Mutations::SequencingSimulations;

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

  if (insert_size==0) {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                ReadSimulator<>::Mode::CREATE, rnd_seed);
  } else {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                insert_size, ReadSimulator<>::Mode::CREATE, rnd_seed);
  }

  simulator.enable_SAM_writing(write_SAM);

  auto mutations_list = forest.get_sample_mutations_list();

  if (FACS) {
    mutations_list = split_by_epigenetic_status(mutations_list, forest);
  }

  auto result = simulator(mutations_list, coverage, purity, "chr_",
                          true, Rcpp::Rcout);

  if (remove_output_path) {
    std::filesystem::remove_all(output_path);
  }

  return get_result_dataframe(result);
}

Rcpp::List simulate_normal_seq(const PhylogeneticForest& forest, const double& coverage, 
                               const int& read_size, const int& insert_size,
                               const std::string& output_dir, const bool& write_SAM,
                               const int& rnd_seed)
{
  using namespace Races::Mutations::SequencingSimulations;

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

  if (insert_size==0) {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                ReadSimulator<>::Mode::CREATE, rnd_seed);
  } else {
    simulator = ReadSimulator<>(output_path, forest.get_reference_path(), read_size,
                                insert_size, ReadSimulator<>::Mode::CREATE, rnd_seed);
  }

  simulator.enable_SAM_writing(write_SAM);

  const auto& germline = forest.get_germline_mutations();

  std::list<Races::Mutations::SampleGenomeMutations> mutations_list;
  mutations_list.emplace_back("Normal sample", germline);
  auto germline_structure_ptr = std::make_shared<Races::Mutations::CellGenomeMutations>(germline.duplicate_structure());
  mutations_list.front().mutations.push_back(germline_structure_ptr);

  auto result = simulator(mutations_list, coverage, 0, "chr_",
                          true, Rcpp::Rcout);

  if (remove_output_path) {
    std::filesystem::remove_all(output_path);
  }

  return get_result_dataframe(result);
}
#include "StateChangePredictor.h"
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
StateChangePredictor::StateChangePredictor(const std::string &predictor_filename,
                                           const cfg::Config &reference_config,
                                           std::set<Element> type_set)
    : type_set_(std::move(type_set)) {
  auto type_set_copy(type_set_);
  type_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(type_set_copy);

  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }
#pragma omp parallel for default(none) shared(reference_config, std::cout)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      auto bond_mapping = GetClusterParametersMappingStateOfBond(reference_config, {i, j});
#pragma omp critical
      {
        site_bond_neighbors_hashmap_[{j, i}] = std::move(bond_mapping);
      }
    }
  }
}
double StateChangePredictor::GetDiffFromAtomIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
double StateChangePredictor::GetDiffFromLatticeIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto [jump_id1, jump_id2] = lattice_id_jump_pair;
  std::pair<size_t, size_t>
      lattice_pair = {std::min(jump_id1, jump_id2), std::max(jump_id1, jump_id2)};
  const auto mapping = site_bond_neighbors_hashmap_.at(lattice_pair);

  const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
  const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);

  int label = 0;
  for (const auto &cluster_vector: mapping) {
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto lattice_id: cluster) {
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id));
        if (lattice_id == lattice_id_jump_pair.first) {
          element_vector_end.push_back(element_second);
          continue;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          element_vector_end.push_back(element_first);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id));
      }
      start_hashmap[cfg::ElementCluster(label, element_vector_start)]++;
      end_hashmap[cfg::ElementCluster(label, element_vector_end)]++;
    }
    label++;
  }

  std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  std::vector<double> de_encode;
  de_encode.reserve(ordered.size());
  const std::vector<double>
      cluster_counter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto count_bond = static_cast<double>(end_hashmap.at(cluster))
        - static_cast<double>(start_hashmap.at(cluster));
    auto total_bond = cluster_counter[static_cast<size_t>(cluster.GetLabel())];
    de_encode.push_back(count_bond / total_bond);
  }
  double dE = 0;
  const size_t cluster_size = base_theta_.size();
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += base_theta_[i] * de_encode[i];
  }
  return dE;
}
} // namespace pred
#include "EnergyPredictorRandomizer.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace pred {
EnergyPredictorRandomizer::EnergyPredictorRandomizer(const std::string &predictor_filename,
                                           const cfg::Config &reference_config,
                                           std::set<Element> type_set)
    : type_set_(std::move(type_set)),
      one_hot_encode_hash_map_(GetOneHotEncodeHashmap(type_set_)),
      mapping_state_(GetClusterParametersMappingState(reference_config)) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  auto type_set_copy(type_set_);
  type_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(type_set_copy);

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }

  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector =
          GetSortedLatticeVectorState(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_state;
      std::transform(sorted_lattice_vector.begin(), sorted_lattice_vector.end(),
                     std::back_inserter(lattice_id_vector_state),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_state_hashmap_[{i, j}] = lattice_id_vector_state;
    }
  }
}
EnergyPredictorRandomizer::~EnergyPredictorRandomizer() = default;

std::pair<double, double> EnergyPredictorRandomizer::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetBarrierAndDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
double EnergyPredictorRandomizer::GetDe(const cfg::Config &config,
                                   const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto &lattice_id_vector = site_bond_cluster_state_hashmap_.at(lattice_id_jump_pair);
  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);

  int label = 0;
  for (const auto &cluster_vector: mapping_state_) {
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto index: cluster) {
        size_t lattice_id = lattice_id_vector[index];
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id));
        if (lattice_id == lattice_id_jump_pair.first) {
          element_vector_end.push_back(migration_element);
          continue;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          element_vector_end.emplace_back(ElementName::X);
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
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto start = static_cast<double>(start_hashmap.at(cluster));
    auto end = static_cast<double>(end_hashmap.at(cluster));
    double total_bond{};
    switch (cluster.GetLabel()) {
      case 0:total_bond = 256;
        break;
      case 1: total_bond = 1536;
        break;
      case 2: total_bond = 768;
        break;
      case 3: total_bond = 3072;
        break;
      case 4: total_bond = 2048;
        break;
      case 5: total_bond = 3072;
        break;
      case 6: total_bond = 6144;
        break;
      case 7: total_bond = 6144;
        break;
      case 8: total_bond = 6144;
        break;
      case 9: total_bond = 6144;
        break;
      case 10: total_bond = 2048;
        break;
    }
    de_encode.push_back((end - start) / total_bond);
  }

  double dE = 0;
  const size_t cluster_size = base_theta_.size();
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += base_theta_[i] * de_encode[i];
  }
  return dE;
}
std::pair<double, double> EnergyPredictorRandomizer::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto dE = GetDe(config, lattice_id_jump_pair);
  return {0, dE};
}
} // namespace pred

#include "EnergyChangePredictorPair.h"
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace pred {
EnergyChangePredictorPair::EnergyChangePredictorPair(const std::string &predictor_filename,
                                                     const cfg::Config &reference_config,
                                                     std::set<Element> element_set)
    : element_set_(std::move(element_set)),
      bond_mapping_state_(GetClusterParametersMappingStatePair(reference_config)) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(element_set_copy);

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      auto base_theta_json = parameters.at("theta");
      base_theta_ = {};
      for (const auto &theta : base_theta_json) {
        base_theta_.emplace_back(theta.get<double>());
      }
      // base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }
#pragma omp parallel for default(none) shared(reference_config)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector = GetSortedLatticeVectorStateOfPair(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_state;
      std::transform(sorted_lattice_vector.begin(), sorted_lattice_vector.end(),
                     std::back_inserter(lattice_id_vector_state),
                     [](const auto &lattice) { return lattice.GetId(); });
#pragma omp critical
      {
        bond_state_hashmap_[{i, j}] = lattice_id_vector_state;
      }
    }
  }
}
EnergyChangePredictorPair::~EnergyChangePredictorPair() = default;
double EnergyChangePredictorPair::GetDeFromAtomIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetDeFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
double EnergyChangePredictorPair::GetDeFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {

  const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
  const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  if (element_first == element_second) {
    return 0.0;
  }
  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);
  const auto &lattice_id_vector = bond_state_hashmap_.at(lattice_id_jump_pair);

#pragma omp parallel for default(none) shared(config, lattice_id_jump_pair, lattice_id_vector, element_first, element_second, start_hashmap, end_hashmap)
  for (size_t label = 0; label < bond_mapping_state_.size(); ++label) {
    const auto &cluster_vector = bond_mapping_state_.at(label);
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto index: cluster) {
        size_t lattice_id = lattice_id_vector[index];
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id));
        if (lattice_id == lattice_id_jump_pair.first) {
          element_vector_end.push_back(element_second);
          continue;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          element_vector_end.emplace_back(element_first);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id));
      }
      auto cluster_start = cfg::ElementCluster(static_cast<int>(label), element_vector_start);
      auto cluster_end = cfg::ElementCluster(static_cast<int>(label), element_vector_end);
#pragma omp critical
      {
        start_hashmap[cluster_start]++;
        end_hashmap[cluster_end]++;
      }
    }
  }

  std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  std::vector<double> de_encode;
  de_encode.reserve(ordered.size());
  static const std::vector<double>
      cluster_counter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto start = static_cast<double>(start_hashmap.at(cluster));
    auto end = static_cast<double>(end_hashmap.at(cluster));
    auto total_bond = cluster_counter[static_cast<size_t>(cluster.GetLabel())];
    de_encode.push_back((end - start) / total_bond);
  }

  double dE = 0;
  const size_t cluster_size = base_theta_.size();
  // not necessary to parallelize this loop
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += base_theta_[i] * de_encode[i];
  }
  return dE;
}
} // pred

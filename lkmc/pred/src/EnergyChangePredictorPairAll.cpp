#include "EnergyChangePredictorPairAll.h"
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace pred {
EnergyChangePredictorPairAll::EnergyChangePredictorPairAll(const std::string &predictor_filename,
                                                           const cfg::Config &reference_config,
                                                           std::set<Element> element_set)
    : element_set_(std::move(element_set)),
      site_mapping_state_(GetClusterParametersMappingStateSite(reference_config)) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(element_set_copy);

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }
#pragma omp parallel for default(none) shared(reference_config)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    auto neighboring_hashset = GetNeighborsLatticeIdSetOfSite(reference_config, i);
    auto sorted_lattice_vector = GetSortedLatticeVectorStateOfSite(reference_config, i);
    std::vector<size_t> lattice_id_vector_state;
    std::transform(sorted_lattice_vector.begin(), sorted_lattice_vector.end(),
                   std::back_inserter(lattice_id_vector_state),
                   [](const auto &lattice) { return lattice.GetId(); });
#pragma omp critical
    {
      site_state_hashmap_[i] = std::move(lattice_id_vector_state);
      neighboring_sites_hashmap_[i] = std::move(neighboring_hashset);
    }
  }
}
EnergyChangePredictorPairAll::~EnergyChangePredictorPairAll() = default;
double EnergyChangePredictorPairAll::GetDeFromAtomIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetDeFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}

double EnergyChangePredictorPairAll::GetDeFromLatticeIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  if (config.GetElementAtLatticeId(lattice_id_jump_pair.first)
      == config.GetElementAtLatticeId(lattice_id_jump_pair.second)) {
    return 0.0;
  }
  const auto &neighbors_of_lhs = neighboring_sites_hashmap_.at(lattice_id_jump_pair.first);
  if (neighbors_of_lhs.find(lattice_id_jump_pair.second) == neighbors_of_lhs.end()) {
    return GetDeFromLatticeIdPairWithoutCoupling(config, lattice_id_jump_pair);
  }
  return GetDeFromLatticeIdPairWithCoupling(config, lattice_id_jump_pair);
}
double EnergyChangePredictorPairAll::GetDeHelper(
    const std::unordered_map<cfg::ElementCluster, size_t,
                             boost::hash<cfg::ElementCluster> > &start_hashmap,
    const std::unordered_map<cfg::ElementCluster, size_t,
                             boost::hash<cfg::ElementCluster> > &end_hashmap,
    const std::map<cfg::ElementCluster, int> &ordered) const {
  std::vector<double> de_encode;
  de_encode.reserve(ordered.size());
  static const std::vector<double>
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
double EnergyChangePredictorPairAll::GetDeFromLatticeIdPairWithCoupling(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
  const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  if (element_first == element_second) {
    return 0.0;
  }
  const auto mapping = GetClusterParametersMappingStatePairOf(config, lattice_id_jump_pair);
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
  return GetDeHelper(start_hashmap, end_hashmap, ordered);
}
double EnergyChangePredictorPairAll::GetDeFromLatticeIdPairWithoutCoupling(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
  const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  if (element_first == element_second) {
    return 0.0;
  }
  auto dE1 = GetDeFromLatticeIdSite(config, lattice_id_jump_pair.first, element_second);
  auto dE2 = GetDeFromLatticeIdSite(config, lattice_id_jump_pair.second, element_first);
  return dE1 + dE2;
}
double EnergyChangePredictorPairAll::GetDeFromLatticeIdSite(const cfg::Config &config,
                                                            size_t lattice_id,
                                                            Element new_element) const {
  const auto old_element = config.GetElementAtLatticeId(lattice_id);
  if (old_element == new_element) {
    return 0.0;
  }

  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);
  const auto &lattice_id_vector = site_state_hashmap_.at(lattice_id);

  int label = 0;
  for (const auto &cluster_vector: site_mapping_state_) {
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto index: cluster) {
        size_t lattice_id_in_cluster = lattice_id_vector[index];
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
        if (lattice_id_in_cluster == lattice_id) {
          element_vector_end.emplace_back(new_element);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
      }
      start_hashmap[cfg::ElementCluster(static_cast<int>(label), element_vector_start)]++;
      end_hashmap[cfg::ElementCluster(static_cast<int>(label), element_vector_end)]++;
    }
    label++;
  }
  std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  return GetDeHelper(start_hashmap, end_hashmap, ordered);
}
} // namespace pred
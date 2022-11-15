#include "EnergyChangePredictorSite.h"
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
EnergyChangePredictorSite::EnergyChangePredictorSite(const std::string &predictor_filename,
                                                     const cfg::Config &reference_config,
                                                     std::set<Element> element_set)
    : element_set_(std::move(element_set)) {
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
    auto site_mapping = GetClusterParametersMappingStateOfLatticeId(reference_config, i);
#pragma omp critical
    {
      site_neighbors_hashmap_[i] = std::move(site_mapping);
    }
  }
}
EnergyChangePredictorSite::~EnergyChangePredictorSite() = default;
double EnergyChangePredictorSite::GetDiffFromAtomId(
    const cfg::Config &config, const size_t atom_id, const Element new_element) const {
  return GetDiffFromLatticeIdPair(config, config.GetLatticeIdFromAtomId(atom_id), new_element);
}
double EnergyChangePredictorSite::GetDiffFromLatticeIdPair(
    const cfg::Config &config, const size_t lattice_id, const Element new_element) const {
  const auto mapping = site_neighbors_hashmap_.at(lattice_id);

  const auto old_element = config.GetElementAtLatticeId(lattice_id);
  if (old_element == new_element) {
    return 0.0;
  }
  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);

  int label = 0;
  for (const auto &cluster_vector: mapping) {
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto lattice_id_in_cluster: cluster) {
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
        if (lattice_id_in_cluster == lattice_id) {
          element_vector_end.push_back(new_element);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
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

// double EnergyChangePredictorSite::GetDiffFromLatticeIdPairSet(
//     const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
//   auto [jump_id1, jump_id2] = lattice_id_jump_pair;
//   std::pair<size_t, size_t>
//       jump_pair = {std::min(jump_id1, jump_id2), std::max(jump_id1, jump_id2)};
//   const auto lattice_id_hashset = site_bond_neighbors_hashmap_.at(jump_pair);
//
//   const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
//   const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
//   auto start_hashmap(initialized_cluster_hashmap_);
//   auto end_hashmap(initialized_cluster_hashmap_);
//
//   for (size_t lattice_id1: lattice_id_hashset) {
//     const Element element1_start = config.GetElementAtLatticeId(lattice_id1);
//     Element element1_end = config.GetElementAtLatticeId(lattice_id1);
//     if (lattice_id1 == lattice_id_jump_pair.first) { element1_end = element_second; }
//     if (lattice_id1 == lattice_id_jump_pair.second) { element1_end = element_first; }
//     start_hashmap[cfg::ElementCluster(0, element1_start)]++;
//     end_hashmap[cfg::ElementCluster(0, element1_end)]++;
//
//     for (size_t lattice_id2: config.GetFirstNeighborsAdjacencyList()[lattice_id1]) {
//       if (lattice_id_hashset.find(lattice_id2) == lattice_id_hashset.end()) { continue; }
//       const Element element2_start = config.GetElementAtLatticeId(lattice_id2);
//       Element element2_end = config.GetElementAtLatticeId(lattice_id2);
//       if (lattice_id2 == lattice_id_jump_pair.first) { element2_end = element_second; }
//       if (lattice_id2 == lattice_id_jump_pair.second) { element2_end = element_first; }
//       start_hashmap[cfg::ElementCluster(1, element1_start, element2_start)]++;
//       end_hashmap[cfg::ElementCluster(1, element1_end, element2_end)]++;
//       for (size_t lattice_id3: config.GetFirstNeighborsAdjacencyList()[lattice_id2]) {
//         if (lattice_id_hashset.find(lattice_id3) == lattice_id_hashset.end()) { continue; }
//         const Element element3_start = config.GetElementAtLatticeId(lattice_id3);
//         Element element3_end = config.GetElementAtLatticeId(lattice_id3);
//         if (lattice_id3 == lattice_id_jump_pair.first) { element3_end = element_second; }
//         if (lattice_id3 == lattice_id_jump_pair.second) { element3_end = element_first; }
//         if (std::find(config.GetFirstNeighborsAdjacencyList()[lattice_id1].begin(),
//                       config.GetFirstNeighborsAdjacencyList()[lattice_id1].end(),
//                       lattice_id3)
//             != config.GetFirstNeighborsAdjacencyList()[lattice_id1].end()) {
//           start_hashmap[cfg::ElementCluster(4, element1_start, element2_start, element3_start)]++;
//           end_hashmap[cfg::ElementCluster(4, element1_end, element2_end, element3_end)]++;
//         } else if (std::find(config.GetSecondNeighborsAdjacencyList()[lattice_id1].begin(),
//                              config.GetSecondNeighborsAdjacencyList()[lattice_id1].end(),
//                              lattice_id3)
//             != config.GetSecondNeighborsAdjacencyList()[lattice_id1].end()) {
//           start_hashmap[cfg::ElementCluster(5, element1_start, element2_start, element3_start)]++;
//           end_hashmap[cfg::ElementCluster(5, element1_end, element2_end, element3_end)]++;
//         } else if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
//                              config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
//                              lattice_id3)
//             != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
//           start_hashmap[cfg::ElementCluster(6, element1_start, element2_start, element3_start)]++;
//           end_hashmap[cfg::ElementCluster(6, element1_end, element2_end, element3_end)]++;
//         }
//       }
//       for (size_t lattice_id3: config.GetSecondNeighborsAdjacencyList()[lattice_id2]) {
//         if (lattice_id_hashset.find(lattice_id3) == lattice_id_hashset.end()) { continue; }
//         const Element element3_start = config.GetElementAtLatticeId(lattice_id3);
//         Element element3_end = config.GetElementAtLatticeId(lattice_id3);
//         if (lattice_id3 == lattice_id_jump_pair.first) { element3_end = element_second; }
//         if (lattice_id3 == lattice_id_jump_pair.second) { element3_end = element_first; }
//         if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
//                       config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
//                       lattice_id3)
//             != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
//           start_hashmap[cfg::ElementCluster(7, element1_start, element2_start, element3_start)]++;
//           end_hashmap[cfg::ElementCluster(7, element1_end, element2_end, element3_end)]++;
//         }
//       }
//     }
//     for (size_t lattice_id2: config.GetSecondNeighborsAdjacencyList()[lattice_id1]) {
//       if (lattice_id_hashset.find(lattice_id2) == lattice_id_hashset.end()) { continue; }
//       const Element element2_start = config.GetElementAtLatticeId(lattice_id2);
//       Element element2_end = config.GetElementAtLatticeId(lattice_id2);
//       if (lattice_id2 == lattice_id_jump_pair.first) { element2_end = element_second; }
//       if (lattice_id2 == lattice_id_jump_pair.second) { element2_end = element_first; }
//       start_hashmap[cfg::ElementCluster(2, element1_start, element2_start)]++;
//       end_hashmap[cfg::ElementCluster(2, element1_end, element2_end)]++;
//     }
//     for (size_t lattice_id2: config.GetThirdNeighborsAdjacencyList()[lattice_id1]) {
//       if (lattice_id_hashset.find(lattice_id2) == lattice_id_hashset.end()) { continue; }
//       const Element element2_start = config.GetElementAtLatticeId(lattice_id2);
//       Element element2_end = config.GetElementAtLatticeId(lattice_id2);
//       if (lattice_id2 == lattice_id_jump_pair.first) { element2_end = element_second; }
//       if (lattice_id2 == lattice_id_jump_pair.second) { element2_end = element_first; }
//       start_hashmap[cfg::ElementCluster(3, element1_start, element2_start)]++;
//       end_hashmap[cfg::ElementCluster(3, element1_end, element2_end)]++;
//     }
//   }
//
//   std::map<cfg::ElementCluster, int>
//       ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
//   std::vector<double> de_encode;
//   de_encode.reserve(ordered.size());
//   const std::vector<double>
//       cluster_counter{256, 3072, 1536, 6144, 12288, 6144, 12288, 6144, 12288, 12288, 12288};
//   for (const auto &cluster_count: ordered) {
//     const auto &cluster = cluster_count.first;
//     auto count_bond = static_cast<double>(end_hashmap.at(cluster))
//         - static_cast<double>(start_hashmap.at(cluster));
//     auto total_bond = cluster_counter[static_cast<size_t>(cluster.GetLabel())];
//     de_encode.push_back(count_bond / total_bond);
//   }
//   double dE = 0;
//   const size_t cluster_size = base_theta_.size();
//   for (size_t i = 0; i < cluster_size; ++i) {
//     dE += base_theta_[i] * de_encode[i];
//   }
//   return dE;
// }
} // namespace pred
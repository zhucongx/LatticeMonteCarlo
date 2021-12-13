#include "ClusterExpansionPredictor.h"

namespace pred {
using Singlet_MM2_t = LatticeClusterMM2<1>;
using Pair_MM2_t = LatticeClusterMM2<2>;
using Triplet_MM2_t = LatticeClusterMM2<3>;
static bool LatticeSortCompare(const cfg::Lattice &lhs,
                               const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double diff_norm = Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
  if (diff_norm < -kEpsilon)
    return true;
  if (diff_norm > kEpsilon)
    return false;
  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  return diff_x < -kEpsilon;
}

template<size_t DataSize>
static bool IsClusterSmallerSymmetrically(const LatticeClusterMM2<DataSize> &lhs,
                                          const LatticeClusterMM2<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    const auto &lhs_lattice = lhs.GetLatticeAt(i);
    const auto &rhs_lattice = rhs.GetLatticeAt(i);
    if (LatticeSortCompare(lhs_lattice, rhs_lattice))
      return true;
    if (LatticeSortCompare(rhs_lattice, lhs_lattice))
      return false;
  }
  // if it reaches here, it means that the clusters are same symmetrically. Returns false.
  return false;
}

// Returns forward and backward sorted lattice lists
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMM2(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  // The number of first-, second-, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = 60;
  auto lattice_id_hashset =
      GetNeighborsLatticeIdSetOfJumpPair(config, lattice_id_jump_pair);
  const auto move_distance = Vector_t{0.5, 0.5, 0.5}
      - GetLatticePairCenter(config, lattice_id_jump_pair);
  std::vector<cfg::Lattice> lattice_list;
  lattice_list.reserve(kNumOfSites);
  for (const auto id: lattice_id_hashset) {
    cfg::Lattice lattice = config.GetLatticeVector()[id];
    // move to center
    auto relative_position = lattice.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);
    lattice.SetRelativePosition(relative_position);
    if (lattice.GetId() == lattice_id_jump_pair.first)
      continue;
    lattice_list.push_back(lattice);
  }
  RotateLatticeVector(lattice_list,
                      GetLatticePairRotationMatrix(config, lattice_id_jump_pair));
  //sort using mmm group point (mm2 if x mirror not applied)
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              if (LatticeSortCompare(lhs, rhs)) return true;
              if (LatticeSortCompare(rhs, lhs)) return false;
              const auto &relative_position_lhs = lhs.GetRelativePosition();
              const auto &relative_position_rhs = rhs.GetRelativePosition();
              const double diff_x =
                  relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
              if (diff_x < -kEpsilon) { return true; }
              if (diff_x > kEpsilon) { return false; }
              const double diff_y =
                  relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
              if (diff_y < -kEpsilon) { return true; }
              if (diff_y > kEpsilon) { return false; }
              const double diff_z =
                  relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
              if (diff_z < -kEpsilon) { return true; }
              if (diff_z > kEpsilon) { return false; }
              return lhs.GetId() < rhs.GetId();
            });
  return lattice_list;
}

template<size_t DataSize>
static void GetAverageParametersMappingFromLatticeClusterVectorHelper(
    std::vector<LatticeClusterMM2<DataSize> > &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  // sort clusters
  std::sort(cluster_vector.begin(), cluster_vector.end(),
            [](const auto &lhs, const auto &rhs) {
              return IsClusterSmallerSymmetrically(lhs, rhs);
            });
  // start to point at Cluster in the first range
  typename std::vector<LatticeClusterMM2<DataSize> >::const_iterator lower_it, upper_it;
  lower_it = cluster_vector.cbegin();
  do {
    upper_it = std::upper_bound(lower_it, cluster_vector.cend(),
                                *lower_it,
                                [](const auto &lhs, const auto &rhs) {
                                  return IsClusterSmallerSymmetrically(lhs, rhs);
                                });
    std::vector<std::vector<size_t> > cluster_index_vector;
    for (auto it = lower_it; it != upper_it; ++it) {
      auto cluster_index = it->GetIndexVector();
      cluster_index_vector.push_back(cluster_index);
    }
    cluster_mapping.push_back(cluster_index_vector);
    // update to next range
    lower_it = upper_it;
  } while (upper_it != cluster_vector.cend());
}

std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMM2(
    const cfg::Config &config) {
  const std::pair<size_t, size_t>
      lattice_id_jump_pair = {0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto lattice_vector = GetSymmetricallySortedLatticeVectorMM2(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  /// singlets
  std::vector<Singlet_MM2_t> singlet_vector;
  for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
    const auto index1 = static_cast<size_t>(std::distance(lattice_vector.begin(), it1));
    cfg::Lattice lattice1(*it1);
    lattice1.SetId(index1);
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::move(singlet_vector), cluster_mapping);
  /// pairs
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > first_pair_set;
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > second_pair_set;
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > third_pair_set;
  for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
    const auto lattice1_index = it1->GetId();
    const auto index1 = static_cast<size_t>(std::distance(lattice_vector.begin(), it1));
    for (const auto &lattice2_index: config.GetFirstNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      const auto index2 = static_cast<size_t>(std::distance(lattice_vector.begin(), it2));
      cfg::Lattice lattice1(*it1), lattice2(*it2);
      lattice1.SetId(index1);
      lattice2.SetId(index2);
      first_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
    }
    for (const auto &lattice2_index: config.GetSecondNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      const auto index2 = static_cast<size_t>(std::distance(lattice_vector.begin(), it2));
      cfg::Lattice lattice1(*it1), lattice2(*it2);
      lattice1.SetId(index1);
      lattice2.SetId(index2);
      second_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
    }
    for (const auto &lattice2_index: config.GetThirdNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      const auto index2 = static_cast<size_t>(std::distance(lattice_vector.begin(), it2));
      cfg::Lattice lattice1(*it1), lattice2(*it2);
      lattice1.SetId(index1);
      lattice2.SetId(index2);
      third_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
    }
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Pair_MM2_t>(first_pair_set.begin(),
                              first_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Pair_MM2_t>(second_pair_set.begin(),
                              second_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Pair_MM2_t>(third_pair_set.begin(),
                              third_pair_set.end()), cluster_mapping);
  // /// triplets
  // std::unordered_set<Triplet_MM2_t, boost::hash<Triplet_MM2_t> > first_first_first_triplets_set;
  // for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
  //   const auto lattice1_index = it1->GetId();
  //   const auto index1 = static_cast<size_t>(std::distance(lattice_vector.begin(), it1));
  //   for (const auto &lattice2_index: config.GetFirstNeighborsAdjacencyList()[lattice1_index]) {
  //     auto it2 = std::find_if(lattice_vector.begin(),
  //                             lattice_vector.end(),
  //                             [lattice2_index](const auto &lattice) {
  //                               return lattice.GetId() == lattice2_index;
  //                             });
  //     if (it2 == lattice_vector.end()) { continue; }
  //     const auto index2 = static_cast<size_t>(std::distance(lattice_vector.begin(), it2));
  //     for (const size_t lattice3_index: config.GetFirstNeighborsAdjacencyList()[lattice2_index]) {
  //       if (lattice3_index == lattice1_index) { continue; }
  //       if (lattice3_index == lattice2_index) { continue; }
  //       if (std::find(config.GetFirstNeighborsAdjacencyList()[lattice1_index].begin(),
  //                     config.GetFirstNeighborsAdjacencyList()[lattice1_index].end(),
  //                     lattice3_index) ==
  //           config.GetFirstNeighborsAdjacencyList()[lattice1_index].end()) { continue; }
  //       auto it3 = std::find_if(lattice_vector.begin(),
  //                               lattice_vector.end(),
  //                               [lattice3_index](const auto &lattice) {
  //                                 return lattice.GetId() == lattice3_index;
  //                               });
  //       if (it3 == lattice_vector.end()) { continue; }
  //       const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
  //       cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
  //       lattice1.SetId(index1);
  //       lattice2.SetId(index2);
  //       lattice3.SetId(index3);
  //       first_first_first_triplets_set.emplace(
  //           std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
  //     }
  //   }
  // }
  // GetAverageParametersMappingFromLatticeClusterVectorHelper(
  //     std::vector<Triplet_MM2_t>(first_first_first_triplets_set.begin(),
  //                                first_first_first_triplets_set.end()), cluster_mapping);

  return cluster_mapping;
}
} // namespace pred
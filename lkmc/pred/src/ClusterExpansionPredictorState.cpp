#include "EnergyPredictor.h"
namespace pred {
using Singlet_State_t = cfg::LatticeCluster<1>;
using Pair_State_t = cfg::LatticeCluster<2>;
using Triplet_State_t = cfg::LatticeCluster<3>;
std::unordered_map<
    cfg::ElementCluster, int, boost::hash<cfg::ElementCluster> > InitializeClusterHashMap(
    const std::set<Element> &type_set) {
  std::unordered_map<cfg::ElementCluster, int,
                     boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap;

  for (const auto &element1: type_set) {
    initialized_cluster_hashmap[cfg::ElementCluster(0, element1)] = 0;
    for (const auto &element2: type_set) {
      if (element2 == ElementName::X) {
        continue;
      }
      if (element1 == ElementName::X && element2.GetString()[0] == 'p') {
        continue;
      }
      for (size_t label = 1; label <= 3; ++label) {
        initialized_cluster_hashmap[cfg::ElementCluster(label, element1, element2)] = 0;
      }
      for (const auto &element3: type_set) {
        if (element3 == ElementName::X || element3.GetString()[0] == 'p') {
          continue;
        }
        for (size_t label = 4; label < 11; ++label) {
          initialized_cluster_hashmap[cfg::ElementCluster(label, element1, element2, element3)] = 0;
        }
      }
    }
  }
  return initialized_cluster_hashmap;
}

static bool LatticeSortCompare(const cfg::Lattice &lhs,
                               const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();
  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  const double diff_y = relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
  const double diff_z = relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
  if (diff_x < -kEpsilon) { return true; }
  if (diff_x > kEpsilon) { return false; }
  if (diff_y < -kEpsilon) { return true; }
  if (diff_y > kEpsilon) { return false; }
  return diff_z < -kEpsilon;
}

// Returns forward and backward sorted lattice lists
std::vector<cfg::Lattice> GetSortedLatticeVectorState(
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
    lattice_list.push_back(lattice);
  }

  RotateLatticeVector(lattice_list,
                      GetLatticePairRotationMatrix(config, lattice_id_jump_pair));
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              return LatticeSortCompare(lhs, rhs);
            });

  return lattice_list;
}

template<size_t DataSize>
static void GetParametersMappingFromLatticeClusterVectorHelper(
    std::vector<cfg::LatticeCluster<DataSize>> &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping
) {
  std::vector<std::vector<size_t> > cluster_index_vector;
  for (const auto &cluster: cluster_vector) {
    auto cluster_index = cluster.GetIndexVector();
    cluster_index_vector.push_back(cluster_index);
  }
  cluster_mapping.push_back(cluster_index_vector);
}
std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingState(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  const auto lattice_vector = GetSortedLatticeVectorState(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  /// singlets
  std::vector<Singlet_State_t> singlet_vector;
  for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{*it1});
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(singlet_vector), cluster_mapping);
  /// pairs
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
  for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
    const auto lattice1_index = it1->GetId();
    for (const auto &lattice2_index: config.GetFirstNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      first_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{*it1, *it2});
    }
    for (const auto &lattice2_index: config.GetSecondNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      second_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{*it1, *it2});
    }
    for (const auto &lattice2_index: config.GetThirdNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      third_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{*it1, *it2});
    }
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(second_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(third_pair_vector), cluster_mapping);
  /// triplets
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;
  std::vector<Triplet_State_t> first_third_third_triplets_vector;
  std::vector<Triplet_State_t> second_third_third_triplets_vector;
  std::vector<Triplet_State_t> third_third_third_triplets_vector;
  for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
    const auto lattice1_index = it1->GetId();
    for (const auto &lattice2_index: config.GetFirstNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      for (const size_t lattice3_index: config.GetFirstNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetFirstNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetFirstNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetFirstNeighborsAdjacencyList()[lattice1_index].end()) {
          first_first_first_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{*it1, *it2, *it3});
        }
        if (std::find(config.GetSecondNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetSecondNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetSecondNeighborsAdjacencyList()[lattice1_index].end()) {
          first_first_second_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{*it1, *it2, *it3});
        }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          first_first_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{*it1, *it2, *it3});
        }
      }
      for (const size_t lattice3_index: config.GetSecondNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          first_second_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{*it1, *it2, *it3});
        }
      }
      for (const size_t lattice3_index: config.GetThirdNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          first_third_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{*it1, *it2, *it3});
        }
      }
    }
    for (const auto &lattice2_index: config.GetSecondNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      for (const size_t lattice3_index: config.GetThirdNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          second_third_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{*it1, *it2, *it3});
        }
      }
    }
    for (const auto &lattice2_index: config.GetThirdNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      for (const size_t lattice3_index: config.GetThirdNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          third_third_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{*it1, *it2, *it3});
        }
      }
    }
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_first_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_second_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_second_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_third_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(second_third_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(third_third_third_triplets_vector), cluster_mapping);
  return cluster_mapping;
}

std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingState(
    const cfg::Config &config) {
  const std::pair<size_t, size_t>
      lattice_id_jump_pair = {0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto lattice_vector = GetSortedLatticeVectorState(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  /// singlets
  std::vector<Singlet_State_t> singlet_vector;
  for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
    const auto index1 = static_cast<size_t>(std::distance(lattice_vector.begin(), it1));
    cfg::Lattice lattice1(*it1);
    lattice1.SetId(index1);
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(singlet_vector), cluster_mapping);
  /// pairs
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
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
      first_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
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
      second_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
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
      third_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
    }
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(second_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(third_pair_vector), cluster_mapping);
  /// triplets
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;
  std::vector<Triplet_State_t> first_third_third_triplets_vector;
  std::vector<Triplet_State_t> second_third_third_triplets_vector;
  std::vector<Triplet_State_t> third_third_third_triplets_vector;
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
      for (const size_t lattice3_index: config.GetFirstNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetFirstNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetFirstNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetFirstNeighborsAdjacencyList()[lattice1_index].end()) {
          const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
          cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
          lattice1.SetId(index1);
          lattice2.SetId(index2);
          lattice3.SetId(index3);
          first_first_first_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
        if (std::find(config.GetSecondNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetSecondNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetSecondNeighborsAdjacencyList()[lattice1_index].end()) {
          const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
          cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
          lattice1.SetId(index1);
          lattice2.SetId(index2);
          lattice3.SetId(index3);
          first_first_second_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
          cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
          lattice1.SetId(index1);
          lattice2.SetId(index2);
          lattice3.SetId(index3);
          first_first_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
      }
      for (const size_t lattice3_index: config.GetSecondNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
          cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
          lattice1.SetId(index1);
          lattice2.SetId(index2);
          lattice3.SetId(index3);
          first_second_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
      }
      for (const size_t lattice3_index: config.GetThirdNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
          cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
          lattice1.SetId(index1);
          lattice2.SetId(index2);
          lattice3.SetId(index3);
          first_third_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
      }
    }
    for (const auto &lattice2_index: config.GetSecondNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      const auto index2 = static_cast<size_t>(std::distance(lattice_vector.begin(), it2));
      for (const size_t lattice3_index: config.GetThirdNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
          cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
          lattice1.SetId(index1);
          lattice2.SetId(index2);
          lattice3.SetId(index3);
          second_third_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
      }
    }
    for (const auto &lattice2_index: config.GetThirdNeighborsAdjacencyList()[lattice1_index]) {
      auto it2 = std::find_if(lattice_vector.begin(),
                              lattice_vector.end(),
                              [lattice2_index](const auto &lattice) {
                                return lattice.GetId() == lattice2_index;
                              });
      if (it2 == lattice_vector.end()) { continue; }
      const auto index2 = static_cast<size_t>(std::distance(lattice_vector.begin(), it2));
      for (const size_t lattice3_index: config.GetThirdNeighborsAdjacencyList()[lattice2_index]) {
        if (lattice3_index == lattice1_index) { continue; }
        if (lattice3_index == lattice2_index) { continue; }
        auto it3 = std::find_if(lattice_vector.begin(),
                                lattice_vector.end(),
                                [lattice3_index](const auto &lattice) {
                                  return lattice.GetId() == lattice3_index;
                                });
        if (it3 == lattice_vector.end()) { continue; }
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_index].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_index].end(),
                      lattice3_index) !=
            config.GetThirdNeighborsAdjacencyList()[lattice1_index].end()) {
          const auto index3 = static_cast<size_t>(std::distance(lattice_vector.begin(), it3));
          cfg::Lattice lattice1(*it1), lattice2(*it2), lattice3(*it3);
          lattice1.SetId(index1);
          lattice2.SetId(index2);
          lattice3.SetId(index3);
          third_third_third_triplets_vector.emplace_back(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
      }
    }
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_first_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_second_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_second_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_third_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(second_third_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(third_third_third_triplets_vector), cluster_mapping);
  return cluster_mapping;
}

std::array<std::vector<double>, 2> GetEncodesFromMapState(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    const std::unordered_map<cfg::ElementCluster,
                             int,
                             boost::hash<cfg::ElementCluster> > &initialized_cluster_hashmap,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto start_hashmap(initialized_cluster_hashmap);
  auto end_hashmap(initialized_cluster_hashmap);
  auto transition_hashmap(initialized_cluster_hashmap);

  size_t label = 0;
  for (const auto &cluster_vector: cluster_mapping) {
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end, element_vector_transition;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      element_vector_transition.reserve(cluster.size());
      for (auto lattice_id: cluster) {
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id));
        if (lattice_id == lattice_id_jump_pair.first) {
          element_vector_end.push_back(migration_element);
          element_vector_transition.push_back(migration_element.GetPseudo());
          continue;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          element_vector_end.emplace_back(ElementName::X);
          element_vector_transition.push_back(migration_element.GetPseudo());
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id));
        element_vector_transition.push_back(config.GetElementAtLatticeId(lattice_id));
      }
      start_hashmap[cfg::ElementCluster(label, element_vector_start)]++;
      end_hashmap[cfg::ElementCluster(label, element_vector_end)]++;
      transition_hashmap[cfg::ElementCluster(label, element_vector_transition)]++;
    }
    label++;
  }
  std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap.begin(), initialized_cluster_hashmap.end());
  std::vector<double> de_encode, e0_encode;
  de_encode.reserve(ordered.size());
  e0_encode.reserve(ordered.size());
  for (const auto &[cluster, count]: ordered) {
    auto start = static_cast<double>(start_hashmap.at(cluster));
    auto end = static_cast<double>(end_hashmap.at(cluster));
    auto transition = static_cast<double>(transition_hashmap.at(cluster));
    double total_bond{};
    switch (cluster.GetLabel()) {
      case 0:total_bond = 256;
        break;
      case 1: total_bond = 3072;
        break;
      case 2: total_bond = 1536;
        break;
      case 3: total_bond = 6144;
        break;
      case 4: total_bond = 12288;
        break;
      case 5: total_bond = 6144;
        break;
      case 6: total_bond = 12288;
        break;
      case 7: total_bond = 6144;
        break;
      case 8: total_bond = 12288;
        break;
      case 9: total_bond = 12288;
        break;
      case 10: total_bond = 12288;
        break;
    }
    de_encode.push_back((end - start) / total_bond);
    e0_encode.push_back((transition - 0.5 * (end + start)) / total_bond);
  }
  return {de_encode, e0_encode};

}

} // namespace pred
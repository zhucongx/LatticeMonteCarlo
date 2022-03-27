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
        for (size_t label = 4; label < 8; ++label) {
          initialized_cluster_hashmap[cfg::ElementCluster(label, element1, element2, element3)] = 0;
        }
      }
    }
  }
  return initialized_cluster_hashmap;
}
static int GetLabel(const std::vector<size_t> &lattice_index_list, const cfg::Config &config) {
  if (lattice_index_list.size() == 1) {
    return 0;
  }
  if (lattice_index_list.size() == 2) {
    return cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[0],
                                                 lattice_index_list[1],
                                                 config);
  }
  if (lattice_index_list.size() == 3) {
    std::vector<int>
        bond_label_list{cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[0],
                                                              lattice_index_list[1],
                                                              config),
                        cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[1],
                                                              lattice_index_list[2],
                                                              config),
                        cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[2],
                                                              lattice_index_list[0],
                                                              config)};
    std::sort(bond_label_list.begin(), bond_label_list.end());
    if (bond_label_list[0] == 1 && bond_label_list[1] == 1 && bond_label_list[2] == 1) {
      return 4;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 1 && bond_label_list[2] == 2) {
      return 5;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 1 && bond_label_list[2] == 3) {
      return 6;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 2 && bond_label_list[2] == 3) {
      return 7;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 3 && bond_label_list[2] == 3) {
      return 8;
    } else if (bond_label_list[0] == 2 && bond_label_list[1] == 3 && bond_label_list[2] == 3) {
      return 9;
    } else if (bond_label_list[0] == 3 && bond_label_list[1] == 3 && bond_label_list[2] == 3) {
      return 10;
    }
  }
  return -1;
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
  std::vector<Singlet_State_t> singlet_vector;
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;
  std::vector<Triplet_State_t> first_third_third_triplets_vector;
  std::vector<Triplet_State_t> second_third_third_triplets_vector;
  std::vector<Triplet_State_t> third_third_third_triplets_vector;

  for (size_t index1 = 0; index1 < lattice_vector.size(); ++index1) {
    cfg::Lattice lattice1(lattice_vector[index1]);
    const size_t id1 = lattice1.GetId();
    lattice1.SetId(index1);
    if (id1 == lattice_id_jump_pair.first || id1 == lattice_id_jump_pair.second) {
      singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    }
    for (size_t index2 = 0; index2 < index1; ++index2) {
      cfg::Lattice lattice2(lattice_vector[index2]);
      const size_t id2 = lattice2.GetId();
      lattice2.SetId(index2);
      if (id1 == lattice_id_jump_pair.first || id1 == lattice_id_jump_pair.second
          || id2 == lattice_id_jump_pair.first || id2 == lattice_id_jump_pair.second) {
        switch (GetLabel({id1, id2}, config)) {
          case 1: first_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 2: second_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 3: third_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          default:continue;
        }
      }
      for (size_t index3 = 0; index3 < index2; ++index3) {
        cfg::Lattice lattice3(lattice_vector[index3]);
        const size_t id3 = lattice3.GetId();
        lattice3.SetId(index3);
        if (id1 == lattice_id_jump_pair.first || id1 == lattice_id_jump_pair.second
            || id2 == lattice_id_jump_pair.first || id2 == lattice_id_jump_pair.second
            || id3 == lattice_id_jump_pair.first || id3 == lattice_id_jump_pair.second) {
          switch (GetLabel({id1, id2, id3}, config)) {
            case 4:
              first_first_first_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 5:
              first_first_second_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 6:
              first_first_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 7:
              first_second_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 8:
              first_third_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 9:
              second_third_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 10:
              third_third_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
          }
        }
      }
    }
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(singlet_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(second_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(third_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_first_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_second_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_second_third_triplets_vector), cluster_mapping);
  // GetParametersMappingFromLatticeClusterVectorHelper(
  //     std::move(first_third_third_triplets_vector), cluster_mapping);
  // GetParametersMappingFromLatticeClusterVectorHelper(
  //     std::move(second_third_third_triplets_vector), cluster_mapping);
  // GetParametersMappingFromLatticeClusterVectorHelper(
  //     std::move(third_third_third_triplets_vector), cluster_mapping);
  return cluster_mapping;
}

} // namespace pred
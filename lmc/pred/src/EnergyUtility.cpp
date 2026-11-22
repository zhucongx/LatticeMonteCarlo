#include "EnergyUtility.h"
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace pred {
std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<Element> &element_set) {
  const size_t type_size = element_set.size();
  std::unordered_map<std::string, std::vector<double> > encode_dict;

  size_t ct1 = 0;
  for (const auto &element: element_set) {
    std::vector<double> element_encode(type_size, 0);
    element_encode[ct1] = 1.0;
    encode_dict[element.GetString()] = element_encode;
    ++ct1;
  }

  const size_t num_pairs = type_size * type_size;
  size_t ct2 = 0;
  for (const auto &element1: element_set) {
    for (const auto &element2: element_set) {
      std::vector<double> element_encode(num_pairs, 0);
      element_encode[ct2] = 1.0;
      encode_dict[element1.GetString() + element2.GetString()] = element_encode;
      ++ct2;
    }
  }

  size_t num_pairs_symmetry = (type_size + 1) * type_size / 2;
  size_t ct3 = 0;
  for (auto it1 = element_set.cbegin(); it1 != element_set.cend(); ++it1) {
    for (auto it2 = it1; it2 != element_set.cend(); ++it2) {
      std::vector<double> element_encode(num_pairs_symmetry, 0);
      element_encode[ct3] = 1.0;
      encode_dict[it1->GetString() + '-' + it2->GetString()] = element_encode;
      ++ct3;
    }
  }

  return encode_dict;
}

std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMMM(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  // The number of the first-nearest neighbors, second-nearest neighbors, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfPair;
  const auto lattice_id_hashset =
      config.GetNeighborsLatticeIdSetOfPair(lattice_id_jump_pair);
  const auto move_distance = Vector_d{0.5, 0.5, 0.5}
      - config.GetLatticePairCenter(lattice_id_jump_pair);
  std::vector<cfg::Lattice> lattice_list;
  lattice_list.reserve(kNumOfSites);
  for (const auto id: lattice_id_hashset) {
    cfg::Lattice lattice = config.GetLatticeVector()[id];
    // move to center
    auto relative_position = lattice.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);
    lattice.SetRelativePosition(relative_position);
    if (lattice.GetId() == lattice_id_jump_pair.first
        || lattice.GetId() == lattice_id_jump_pair.second)
      continue;
    lattice_list.push_back(lattice);
  }
  RotateLatticeVector(lattice_list,
                      config.GetLatticePairRotationMatrix(lattice_id_jump_pair));
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              return PositionCompareMMM(lhs, rhs);
            });
  return lattice_list;
}
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMM2(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  // The number of first-nearest neighbors, second-nearest neighbors, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfPair;
  const auto lattice_id_hashset =
      config.GetNeighborsLatticeIdSetOfPair(lattice_id_jump_pair);
  const auto move_distance = Vector_d{0.5, 0.5, 0.5}
      - config.GetLatticePairCenter(lattice_id_jump_pair);
  std::vector<cfg::Lattice> lattice_list;
  lattice_list.reserve(kNumOfSites);
  for (const auto id: lattice_id_hashset) {
    cfg::Lattice lattice = config.GetLatticeVector()[id];
    // move to center
    auto relative_position = lattice.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);
    lattice.SetRelativePosition(relative_position);
    if (lattice.GetId() == lattice_id_jump_pair.first
        || lattice.GetId() == lattice_id_jump_pair.second)
      continue;
    lattice_list.push_back(lattice);
  }
  RotateLatticeVector(lattice_list,
                      config.GetLatticePairRotationMatrix(lattice_id_jump_pair));
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              return PositionCompareMM2(lhs, rhs);
            });
  return lattice_list;
}

template<size_t DataSize>
static void GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
    std::vector<cfg::LatticeClusterMMM<DataSize> > &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  // sort clusters
  std::sort(cluster_vector.begin(), cluster_vector.end(),
            [](const auto &lhs, const auto &rhs) -> bool {
              return cfg::IsClusterSmallerSymmetricallyMMM(lhs, rhs);
            });
  // start to point at Cluster in the first range
  typename std::vector<cfg::LatticeClusterMMM<DataSize> >::const_iterator lower_it, upper_it;
  lower_it = cluster_vector.cbegin();
  do {
    upper_it = std::upper_bound(lower_it, cluster_vector.cend(),
                                *lower_it,
                                [](const auto &lhs, const auto &rhs) {
                                  return cfg::IsClusterSmallerSymmetricallyMMM(lhs, rhs);
                                });
    std::vector<std::vector<size_t> > cluster_index_vector;
    for (auto it = lower_it; it != upper_it; ++it) {
      auto cluster_index = it->GetIndexVector();
      if (it->IsSymmetryLabel()) {
        cluster_index.insert(cluster_index.begin(), SIZE_MAX);
      }
      cluster_index_vector.push_back(cluster_index);
    }
    cluster_mapping.push_back(cluster_index_vector);
    // update to next range
    lower_it = upper_it;
  } while (upper_it != cluster_vector.cend());
}
template<size_t DataSize>
static void GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
    std::vector<cfg::LatticeClusterMM2<DataSize> > &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  // sort clusters
  std::sort(cluster_vector.begin(), cluster_vector.end(),
            [](const auto &lhs, const auto &rhs) {
              return cfg::IsClusterSmallerSymmetricallyMM2(lhs, rhs);
            });
  // start to point at Cluster in the first range
  typename std::vector<cfg::LatticeClusterMM2<DataSize> >::const_iterator lower_it, upper_it;
  lower_it = cluster_vector.cbegin();
  do {
    upper_it = std::upper_bound(lower_it, cluster_vector.cend(),
                                *lower_it,
                                [](const auto &lhs, const auto &rhs) {
                                  return cfg::IsClusterSmallerSymmetricallyMM2(lhs, rhs);
                                });
    std::vector<std::vector<size_t> > cluster_index_vector;
    for (auto it = lower_it; it != upper_it; ++it) {
      auto cluster_index = it->GetIndexVector();
      if (it->IsSymmetryLabel()) {
        cluster_index.insert(cluster_index.begin(), SIZE_MAX);
      }
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
  // singlets
  std::vector<Singlet_MM2_t> singlet_vector;
  // pairs
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > first_pair_set;
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > second_pair_set;
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > third_pair_set;
  for (size_t index1 = 0; index1 < lattice_vector.size(); ++index1) {
    cfg::Lattice lattice1(lattice_vector[index1]);
    const size_t id1 = lattice1.GetId();
    lattice1.SetId(index1);
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    for (size_t index2 = 0; index2 < index1; ++index2) {
      cfg::Lattice lattice2(lattice_vector[index2]);
      const size_t id2 = lattice2.GetId();
      lattice2.SetId(index2);
      switch (config.FindDistanceLabelBetweenLattice(id1, id2)) {
        case 1: first_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 2: second_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 3: third_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        default: break; // Added to resolve -Wswitch-default
      }
    }
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::move(singlet_vector), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::vector<Pair_MM2_t>(first_pair_set.begin(),
                              first_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::vector<Pair_MM2_t>(second_pair_set.begin(),
                              second_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::vector<Pair_MM2_t>(third_pair_set.begin(),
                              third_pair_set.end()), cluster_mapping);

  return cluster_mapping;
}
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMMM(
    const cfg::Config &config) {
  const std::pair<size_t, size_t>
      lattice_id_jump_pair = {0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto lattice_vector = GetSymmetricallySortedLatticeVectorMMM(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  // singlets
  std::vector<Singlet_MMM_t> singlet_vector;
  // pairs
  std::unordered_set<Pair_MMM_t, boost::hash<Pair_MMM_t> > first_pair_set;
  std::unordered_set<Pair_MMM_t, boost::hash<Pair_MMM_t> > second_pair_set;
  std::unordered_set<Pair_MMM_t, boost::hash<Pair_MMM_t> > third_pair_set;
  for (size_t index1 = 0; index1 < lattice_vector.size(); ++index1) {
    cfg::Lattice lattice1(lattice_vector[index1]);
    const size_t id1 = lattice1.GetId();
    lattice1.SetId(index1);
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    for (size_t index2 = 0; index2 < index1; ++index2) {
      cfg::Lattice lattice2(lattice_vector[index2]);
      const size_t id2 = lattice2.GetId();
      lattice2.SetId(index2);
      switch (config.FindDistanceLabelBetweenLattice(id1, id2)) {
        case 1: first_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 2: second_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 3: third_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        default: break; // Added to resolve -Wswitch-default
      }
    }
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::move(singlet_vector), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::vector<Pair_MMM_t>(first_pair_set.begin(),
                              first_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::vector<Pair_MMM_t>(second_pair_set.begin(),
                              second_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::vector<Pair_MMM_t>(third_pair_set.begin(),
                              third_pair_set.end()), cluster_mapping);
  return cluster_mapping;
}

std::vector<cfg::Lattice> GetSortedLatticeVectorStateOfPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_pair) {
  // The number of first-nearest neighbors, second-nearest neighbors, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfPair;
  const auto lattice_id_hashset = config.GetNeighborsLatticeIdSetOfPair(lattice_id_pair);
  const auto move_distance = Vector_d{0.5, 0.5, 0.5} - config.GetLatticePairCenter(lattice_id_pair);
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
                      config.GetLatticePairRotationMatrix(lattice_id_pair));
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              return PositionCompareState(lhs, rhs);
            });

  return lattice_list;
}
std::vector<cfg::Lattice> GetSortedLatticeVectorStateOfSite(
    const cfg::Config &config, const size_t lattice_id) {
  // The number of the first-nearest neighbors, second-nearest neighbors, and third-nearest neighbors of the lattice sites
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfSite;
  const auto lattice_id_hashset = config.GetNeighborsLatticeIdSetOfSite(lattice_id);
  const auto move_distance =
      Vector_d{0.5, 0.5, 0.5} - config.GetLatticeVector()[lattice_id].GetRelativePosition();
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

  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              return PositionCompareState(lhs, rhs);
            });

  return lattice_list;
}
std::unordered_map<
    cfg::ElementCluster, size_t, boost::hash<cfg::ElementCluster> > InitializeClusterHashMap(
    const std::set<Element> &element_set) {
  std::unordered_map<cfg::ElementCluster, size_t,
                     boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap;

  for (const auto &element1: element_set) {
    initialized_cluster_hashmap[cfg::ElementCluster(0, element1)] = 0;
    for (const auto &element2: element_set) {
      if (element2 == ElementName::X) {
        continue;
      }
      if (element1 == ElementName::X && element2.GetString()[0] == 'p') {
        continue;
      }
      for (int label = 1; label <= 3; ++label) {
        initialized_cluster_hashmap[cfg::ElementCluster(label, element1, element2)] = 0;
      }
      for (const auto &element3: element_set) {
        if (element3 == ElementName::X || element3.GetString()[0] == 'p') {
          continue;
        }
        for (int label = 4; label < 8; ++label) {
          initialized_cluster_hashmap[cfg::ElementCluster(label, element1, element2, element3)] = 0;
        }
      }
    }
  }
  return initialized_cluster_hashmap;
}

int GetLabel(const std::vector<size_t> &lattice_index_list, const cfg::Config &config) {
  if (lattice_index_list.size() == 1) {
    return 0;
  }
  if (lattice_index_list.size() == 2) {
    return config.FindDistanceLabelBetweenLattice(lattice_index_list[0],
                                                  lattice_index_list[1]);
  }
  if (lattice_index_list.size() == 3) {
    std::vector<int>
        bond_label_list{config.FindDistanceLabelBetweenLattice(lattice_index_list[0],
                                                               lattice_index_list[1]),
                        config.FindDistanceLabelBetweenLattice(lattice_index_list[1],
                                                               lattice_index_list[2]),
                        config.FindDistanceLabelBetweenLattice(lattice_index_list[2],
                                                               lattice_index_list[0])};
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

template<size_t DataSize>
static void GetParametersMappingFromLatticeClusterVectorHelper(
    std::vector<cfg::LatticeCluster<DataSize>> &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  std::vector<std::vector<size_t> > cluster_index_vector;
  for (const auto &cluster: cluster_vector) {
    auto cluster_index = cluster.GetIndexVector();
    cluster_index_vector.push_back(cluster_index);
  }
  cluster_mapping.push_back(cluster_index_vector);
}

std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingStatePair(
    const cfg::Config &config) {
  const std::pair<size_t, size_t>
      lattice_id_jump_pair = {0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto lattice_vector = GetSortedLatticeVectorStateOfPair(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  std::vector<Singlet_State_t> singlet_vector;
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;
  // std::vector<Triplet_State_t> first_third_third_triplets_vector;
  // std::vector<Triplet_State_t> second_third_third_triplets_vector;
  // std::vector<Triplet_State_t> third_third_third_triplets_vector;

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
              // first_third_third_triplets_vector.emplace_back(
              //     std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 9:
              // second_third_third_triplets_vector.emplace_back(
              //     std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 10:
              // third_third_third_triplets_vector.emplace_back(
              //     std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            default:continue;
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
std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingStateSite(
    const cfg::Config &config) {
  const size_t lattice_id = 0;
  const auto lattice_vector = GetSortedLatticeVectorStateOfSite(config, lattice_id);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  std::vector<Singlet_State_t> singlet_vector;
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;

  for (size_t index1 = 0; index1 < lattice_vector.size(); ++index1) {
    cfg::Lattice lattice1(lattice_vector[index1]);
    const size_t id1 = lattice1.GetId();
    lattice1.SetId(index1);
    if (id1 == lattice_id) {
      singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    }
    for (size_t index2 = 0; index2 < index1; ++index2) {
      cfg::Lattice lattice2(lattice_vector[index2]);
      const size_t id2 = lattice2.GetId();
      lattice2.SetId(index2);
      if (id1 == lattice_id || id2 == lattice_id) {
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
        if (id1 == lattice_id || id2 == lattice_id || id3 == lattice_id) {
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
            default:continue;
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
  return cluster_mapping;
}
std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingStatePairOf(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_pair) {
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  std::vector<Singlet_State_t> singlet_vector;
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;

  auto lattice_id_hashset = config.GetNeighborsLatticeIdSetOfPair(lattice_id_pair);

  for (auto it1 = lattice_id_hashset.begin(); it1 != lattice_id_hashset.end(); ++it1) {
    auto index1 = *it1;
    const auto &lattice1 = config.GetLatticeVector()[index1];
    if (index1 == lattice_id_pair.first || index1 == lattice_id_pair.second) {
      singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    }
    for (auto it2 = lattice_id_hashset.begin(); it2 != it1; ++it2) {
      auto index2 = *it2;
      const auto &lattice2 = config.GetLatticeVector()[index2];
      if (index1 == lattice_id_pair.first || index1 == lattice_id_pair.second
          || index2 == lattice_id_pair.first || index2 == lattice_id_pair.second) {
        switch (GetLabel({index1, index2}, config)) {
          case 1: first_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 2: second_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 3: third_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          default:continue;
        }
      }
      for (auto it3 = lattice_id_hashset.begin(); it3 != it2; ++it3) {
        auto index3 = *it3;
        const auto &lattice3 = config.GetLatticeVector()[index3];
        if (index1 == lattice_id_pair.first || index1 == lattice_id_pair.second
            || index2 == lattice_id_pair.first || index2 == lattice_id_pair.second
            || index3 == lattice_id_pair.first || index3 == lattice_id_pair.second) {
          switch (GetLabel({index1, index2, index3}, config)) {
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
            default:continue;
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
  return cluster_mapping;
}
std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingStateSiteOf(
    const cfg::Config &config, size_t lattice_id) {
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  std::vector<Singlet_State_t> singlet_vector;
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;

  auto lattice_id_hashset = config.GetNeighborsLatticeIdSetOfSite(lattice_id);

  for (auto it1 = lattice_id_hashset.begin(); it1 != lattice_id_hashset.end(); ++it1) {
    auto index1 = *it1;
    const auto &lattice1 = config.GetLatticeVector()[index1];
    if (index1 == lattice_id) {
      singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    }
    for (auto it2 = lattice_id_hashset.begin(); it2 != it1; ++it2) {
      auto index2 = *it2;
      const auto &lattice2 = config.GetLatticeVector()[index2];
      if (index1 == lattice_id || index2 == lattice_id) {
        switch (GetLabel({index1, index2}, config)) {
          case 1: first_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 2: second_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 3: third_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          default:continue;
        }
      }
      for (auto it3 = lattice_id_hashset.begin(); it3 != it2; ++it3) {
        auto index3 = *it3;
        const auto &lattice3 = config.GetLatticeVector()[index3];
        if (index1 == lattice_id || index2 == lattice_id || index3 == lattice_id) {
          switch (GetLabel({index1, index2, index3}, config)) {
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
            default:continue;
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
  return cluster_mapping;
}
// Number of one-hot slots taken by one group of symmetry-equivalent clusters.
static size_t GetOneHotListLength(const size_t num_of_elements,
                                  const std::vector<std::vector<size_t>> &cluster_vector) {
  if (cluster_vector[0][0] == SIZE_MAX) {
    return num_of_elements * (num_of_elements + 1) / 2;
  }
  size_t list_length = 1;
  for (size_t i = 0; i < cluster_vector[0].size(); ++i) {
    list_length *= num_of_elements;
  }
  return list_length;
}

size_t GetOneHotEncodeSize(const size_t num_of_elements,
                           const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping) {
  size_t encode_size = 0;
  for (const auto &cluster_vector: cluster_mapping) {
    encode_size += GetOneHotListLength(num_of_elements, cluster_vector);
  }
  return encode_size;
}

void CheckModelSize(const std::string &predictor_filename,
                    const std::string &name,
                    const size_t actual_size,
                    const size_t expected_size) {
  if (actual_size != expected_size) {
    throw std::runtime_error(predictor_filename + ": " + name + " has " + std::to_string(actual_size) +
                             " entries but " + std::to_string(expected_size) + " are expected. " +
                             "The coefficient file does not match element_set or is inconsistent.");
  }
}

void GetOneHotParametersFromMap(
    const std::vector<Element> &encode,
    const std::unordered_map<std::string, std::vector<double>> &one_hot_encode_hashmap,
    const size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping,
    std::vector<double> &out_encode) {
  out_encode.clear();
  // out_encode.reserve(2452);
  // Todo 711 for mmm ternary, 1401 for mmm2 ternary, 1240 for mmm quaternary, 2452 for mmm2 quaternary
  auto &sum_of_list = GetThreadLocalSumOfListBuffer();
  for (const auto &cluster_vector: cluster_mapping) {
    const size_t list_length = GetOneHotListLength(num_of_elements, cluster_vector);
    sum_of_list.assign(list_length, 0.0);

    for (const auto &cluster: cluster_vector) {
      std::string cluster_type;
      if (cluster[0] == SIZE_MAX) {
        std::vector<std::string> string_vector;
        std::transform(std::next(cluster.begin()), cluster.end(),
                       std::back_inserter(string_vector),
                       [&encode](const auto &index) { return encode[index].GetString(); });
        std::sort(string_vector.begin(), string_vector.end());
        cluster_type = string_vector.empty() ? "" :
                       std::accumulate(
                           std::next(string_vector.begin()), string_vector.end(),
                           *string_vector.begin(),
                           [](auto &&a, auto &&b) -> auto & {
                             a += '-';
                             a += b;
                             return a;
                           });
      } else {
        for (const auto index: cluster) {
          cluster_type += encode[index].GetString();
        }
      }
      const auto &cluster_one_hot_encode = one_hot_encode_hashmap.at(cluster_type);
      std::transform(sum_of_list.begin(), sum_of_list.end(),
                     cluster_one_hot_encode.begin(),
                     sum_of_list.begin(),
                     std::plus<>());
    }
    auto cluster_vector_size = static_cast<double>( cluster_vector.size());
    std::for_each(sum_of_list.begin(), sum_of_list.end(),
                  [cluster_vector_size](auto &n) { n /= cluster_vector_size; });

    std::move(sum_of_list.begin(), sum_of_list.end(), std::back_inserter(out_encode));
  }
}

FlatClusterMapping FlattenClusterMapping(const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping) {
  FlatClusterMapping flat_mapping;
  for (size_t label = 0; label < cluster_mapping.size(); ++label) {
    const auto &cluster_vector = cluster_mapping[label];
    if (cluster_vector.empty()) {
      continue;
    }
    const size_t cluster_size = cluster_vector.front().size();
    if (cluster_size == 0 || cluster_size > FlatClusterMapping::kMaxClusterSize) {
      throw std::invalid_argument("A cluster must have between 1 and 3 sites");
    }
    flat_mapping.groups.push_back({label, cluster_size, cluster_vector.size(), flat_mapping.site_indices.size()});
    for (const auto &cluster: cluster_vector) {
      if (cluster.size() != cluster_size) {
        throw std::invalid_argument("Clusters with the same label must have the same number of sites");
      }
      for (const auto site_index: cluster) {
        if (site_index > std::numeric_limits<uint8_t>::max()) {
          throw std::invalid_argument("Site index of a cluster is too large for the flat cluster mapping");
        }
        flat_mapping.site_indices.push_back(static_cast<uint8_t>(site_index));
      }
    }
  }
  return flat_mapping;
}

OneHotEncoder::OneHotEncoder(const std::set<Element> &element_set,
                             const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping)
    : num_elements_(element_set.size()) {
  // The position of an element in the set is its position in every block, as in GetOneHotEncodeHashmap().
  element_indices_.fill(kNoElement);
  uint8_t element_index = 0;
  for (const auto element: element_set) {
    element_indices_[static_cast<uint8_t>(hash_value(element))] = element_index++;
  }
  // Pairs of equivalent sites are listed as (0, 0), (0, 1), ..., (1, 1), (1, 2), ...
  symmetric_row_offsets_.resize(num_elements_);
  size_t num_symmetric_pairs = 0;
  for (size_t i = 0; i < num_elements_; ++i) {
    symmetric_row_offsets_[i] = num_symmetric_pairs - i;
    num_symmetric_pairs += num_elements_ - i;
  }

  for (const auto &cluster_vector: cluster_mapping) {
    if (cluster_vector.empty()) {
      throw std::invalid_argument("A group of the one-hot cluster mapping has no cluster");
    }
    // A cluster of equivalent sites starts with SIZE_MAX, see GetAverageClusterParametersMappingMMM().
    const bool is_symmetric = cluster_vector.front().front() == SIZE_MAX;
    const size_t cluster_size = cluster_vector.front().size() - (is_symmetric ? 1 : 0);
    if (cluster_size == 0 || cluster_size > 2 || (is_symmetric && cluster_size != 2)) {
      throw std::invalid_argument("The one-hot encoding is defined for single sites and pairs of sites only");
    }
    Group group;
    group.cluster_size = cluster_size;
    group.is_symmetric = is_symmetric;
    group.num_clusters = cluster_vector.size();
    group.site_offset = site_indices_.size();
    group.encode_offset = encode_size_;
    group.encode_length = is_symmetric ? num_symmetric_pairs : (cluster_size == 1 ? num_elements_
                                                                                  : num_elements_ * num_elements_);
    for (const auto &cluster: cluster_vector) {
      const bool is_cluster_symmetric = cluster.front() == SIZE_MAX;
      if (is_cluster_symmetric != is_symmetric || cluster.size() != cluster_size + (is_symmetric ? 1 : 0)) {
        throw std::invalid_argument("The clusters of a group of the one-hot cluster mapping differ in kind");
      }
      for (size_t i = is_symmetric ? 1 : 0; i < cluster.size(); ++i) {
        if (cluster[i] > std::numeric_limits<uint8_t>::max()) {
          throw std::invalid_argument("Site index of a cluster is too large for the one-hot encoding");
        }
        site_indices_.push_back(static_cast<uint8_t>(cluster[i]));
      }
    }
    encode_size_ += group.encode_length;
    groups_.push_back(group);
  }
}

void OneHotEncoder::ThrowElementNotFound(const Element element) {
  throw std::out_of_range("Element " + element.GetString() + " is not part of the one-hot encoding");
}

void OneHotEncoder::Encode(const uint8_t *element_indices, double *out_encode) const {
  for (const auto &group: groups_) {
    double *block = out_encode + group.encode_offset;
    std::fill_n(block, group.encode_length, 0.0);
    const uint8_t *sites = site_indices_.data() + group.site_offset;
    if (group.cluster_size == 1) {
      for (size_t n = 0; n < group.num_clusters; ++n) {
        block[element_indices[sites[n]]] += 1.0;
      }
    } else if (group.is_symmetric) {
      for (size_t n = 0; n < group.num_clusters; ++n, sites += 2) {
        const size_t lower = std::min(element_indices[sites[0]], element_indices[sites[1]]);
        const size_t upper = std::max(element_indices[sites[0]], element_indices[sites[1]]);
        block[symmetric_row_offsets_[lower] + upper] += 1.0;
      }
    } else {
      for (size_t n = 0; n < group.num_clusters; ++n, sites += 2) {
        block[element_indices[sites[0]] * num_elements_ + element_indices[sites[1]]] += 1.0;
      }
    }
    const auto num_clusters = static_cast<double>(group.num_clusters);
    for (size_t i = 0; i < group.encode_length; ++i) {
      block[i] /= num_clusters;
    }
  }
}

std::vector<double> OneHotEncoder::ScaleWeights(const Eigen::VectorXd &weights) const {
  if (static_cast<size_t>(weights.size()) != encode_size_) {
    throw std::invalid_argument("The weights do not have the size of the one-hot encoding");
  }
  std::vector<double> scaled_weights(encode_size_);
  for (const auto &group: groups_) {
    const auto num_clusters = static_cast<double>(group.num_clusters);
    for (size_t i = group.encode_offset; i < group.encode_offset + group.encode_length; ++i) {
      scaled_weights[i] = weights[static_cast<Eigen::Index>(i)] / num_clusters;
    }
  }
  return scaled_weights;
}

double OneHotEncoder::Dot(const uint8_t *element_indices, const double *scaled_weights) const {
  double sum = 0.0;
  for (const auto &group: groups_) {
    const double *block = scaled_weights + group.encode_offset;
    const uint8_t *sites = site_indices_.data() + group.site_offset;
    double group_sum = 0.0;
    if (group.cluster_size == 1) {
      for (size_t n = 0; n < group.num_clusters; ++n) {
        group_sum += block[element_indices[sites[n]]];
      }
    } else if (group.is_symmetric) {
      for (size_t n = 0; n < group.num_clusters; ++n, sites += 2) {
        const size_t lower = std::min(element_indices[sites[0]], element_indices[sites[1]]);
        const size_t upper = std::max(element_indices[sites[0]], element_indices[sites[1]]);
        group_sum += block[symmetric_row_offsets_[lower] + upper];
      }
    } else {
      for (size_t n = 0; n < group.num_clusters; ++n, sites += 2) {
        group_sum += block[element_indices[sites[0]] * num_elements_ + element_indices[sites[1]]];
      }
    }
    sum += group_sum;
  }
  return sum;
}

ClusterIndexer::ClusterIndexer(const std::map<cfg::ElementCluster, int> &ordered_clusters,
                               std::vector<double> cluster_total_bonds)
    : cluster_total_bonds_(std::move(cluster_total_bonds)) {
  cluster_index_to_cluster_.reserve(ordered_clusters.size());
  for (const auto &entry : ordered_clusters) {
    const auto &cluster = entry.first;
    const size_t idx = cluster_index_to_cluster_.size();
    cluster_index_to_cluster_.push_back(cluster);
    cluster_to_index_[cluster] = idx;
  }
  if (cluster_index_to_cluster_.size() != cluster_total_bonds_.size()) {
    throw std::invalid_argument("cluster_total_bonds size mismatch in ClusterIndexer");
  }
  BuildLookup();
}

void ClusterIndexer::BuildLookup() {
  // Every element that is part of a cluster gets a code of its own, all other elements share the code after them.
  constexpr uint8_t kNoCode = std::numeric_limits<uint8_t>::max();
  element_codes_.fill(kNoCode);
  num_codes_ = 0;
  size_t num_labels = 0;
  for (const auto &cluster: cluster_index_to_cluster_) {
    if (cluster.GetLabel() < 0 || cluster.GetSize() == 0 || cluster.GetSize() > FlatClusterMapping::kMaxClusterSize) {
      throw std::invalid_argument("A cluster needs a non-negative label and between 1 and 3 elements");
    }
    num_labels = std::max(num_labels, static_cast<size_t>(cluster.GetLabel()) + 1);
    for (const auto element: cluster.GetElementVector()) {
      auto &code = element_codes_[static_cast<uint8_t>(hash_value(element))];
      if (code == kNoCode) {
        code = static_cast<uint8_t>(num_codes_++);
      }
    }
  }
  for (auto &code: element_codes_) {
    if (code == kNoCode) {
      code = static_cast<uint8_t>(num_codes_);
    }
  }
  ++num_codes_;
  codes_per_label_ = num_codes_ * num_codes_ * num_codes_;

  // Enter every cluster under each order of its elements, so that a lookup does not have to sort them.
  label_sizes_.assign(num_labels, std::numeric_limits<size_t>::max());
  lookup_.assign(num_labels * codes_per_label_, kNotFound);
  for (size_t index = 0; index < cluster_index_to_cluster_.size(); ++index) {
    const auto &cluster = cluster_index_to_cluster_[index];
    const auto label = static_cast<size_t>(cluster.GetLabel());
    const auto size = cluster.GetSize();
    if (label_sizes_[label] != std::numeric_limits<size_t>::max() && label_sizes_[label] != size) {
      throw std::invalid_argument("Clusters with the same label must have the same number of elements");
    }
    label_sizes_[label] = size;
    std::array<uint8_t, FlatClusterMapping::kMaxClusterSize> codes{};
    for (size_t i = 0; i < size; ++i) {
      codes[i] = GetElementCode(cluster.GetElementVector()[i]);
    }
    const auto codes_end = codes.begin() + static_cast<std::ptrdiff_t>(size);
    std::sort(codes.begin(), codes_end);
    do {
      size_t combined_code = 0;
      for (size_t i = 0; i < size; ++i) {
        combined_code = combined_code * num_codes_ + codes[i];
      }
      lookup_[label * codes_per_label_ + combined_code] = static_cast<uint32_t>(index);
    } while (std::next_permutation(codes.begin(), codes_end));
  }
}

void ClusterIndexer::ThrowClusterNotFound() {
  throw std::out_of_range("Cluster not found in ClusterIndexer");
}

size_t ClusterIndexer::GetIndex(const cfg::ElementCluster &cluster) const {
  const auto it = cluster_to_index_.find(cluster);
  if (it == cluster_to_index_.end()) {
    ThrowClusterNotFound();
  }
  return it->second;
}
} // pred

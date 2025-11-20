#include "EnergyUtility.h"
namespace pred {
std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<Element> &element_set) {
  size_t type_size = element_set.size();
  std::unordered_map<std::string, std::vector<double> > encode_dict;

  size_t ct1 = 0;
  for (const auto &element: element_set) {
    std::vector<double> element_encode(type_size, 0);
    element_encode[ct1] = 1.0;
    encode_dict[element.GetString()] = element_encode;
    ++ct1;
  }

  size_t num_pairs = type_size * type_size;
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
  // The number of first-, second-, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfPair;
  auto lattice_id_hashset =
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
  // The number of first-, second-, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfPair;
  auto lattice_id_hashset =
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
  // The number of first-, second-, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfPair;
  auto lattice_id_hashset =
      config.GetNeighborsLatticeIdSetOfPair(lattice_id_pair);
  const auto move_distance = Vector_d{0.5, 0.5, 0.5}
      - config.GetLatticePairCenter(lattice_id_pair);
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
  // The number of first-, second-, and third-nearest neighbors of the lattice sites
  constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfSite;
  auto lattice_id_hashset = config.GetNeighborsLatticeIdSetOfSite(lattice_id);
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
std::vector<double> GetOneHotParametersFromMap(
    const std::vector<Element> &encode,
    const std::unordered_map<std::string, std::vector<double> > &one_hot_encode_hashmap,
    size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {

  std::vector<double> res_encode;
  res_encode.reserve(1401); // Todo 711 for mmm and 1401 for mmm2
  for (const auto &cluster_vector: cluster_mapping) {
    size_t list_length;
    if (cluster_vector[0][0] == SIZE_MAX) {
      list_length = static_cast<size_t>(num_of_elements * (num_of_elements + 1) / 2);
    } else {
      list_length = static_cast<size_t>(std::pow(num_of_elements, cluster_vector[0].size()));
    }
    std::vector<double> sum_of_list(list_length, 0);

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
        for (auto index: cluster) {
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
    std::for_each(sum_of_list.begin(),
                  sum_of_list.end(),
                  [cluster_vector_size](auto &n) { n /= cluster_vector_size; });

    std::move(sum_of_list.begin(), sum_of_list.end(), std::back_inserter(res_encode));
  }
  return res_encode;
}
} // pred

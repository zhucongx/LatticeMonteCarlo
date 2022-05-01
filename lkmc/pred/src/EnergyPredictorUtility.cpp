#include "EnergyPredictorUtility.hpp"
namespace pred {
std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<Element> &type_set) {
  size_t type_size = type_set.size();
  std::unordered_map<std::string, std::vector<double> > encode_dict;

  size_t ct1 = 0;
  for (const auto &element: type_set) {
    std::vector<double> element_encode(type_size, 0);
    element_encode[ct1] = 1.0;
    encode_dict[element.GetString()] = element_encode;
    ++ct1;
  }

  size_t num_pairs = type_size * type_size;
  size_t ct2 = 0;
  for (const auto &element1: type_set) {
    for (const auto &element2: type_set) {
      std::vector<double> element_encode(num_pairs, 0);
      element_encode[ct2] = 1.0;
      encode_dict[element1.GetString() + element2.GetString()] = element_encode;
      ++ct2;
    }
  }

  size_t num_pairs_symmetry = (type_size + 1) * type_size / 2;
  size_t ct3 = 0;
  for (auto it1 = type_set.cbegin(); it1 != type_set.cend(); ++it1) {
    for (auto it2 = it1; it2 != type_set.cend(); ++it2) {
      std::vector<double> element_encode(num_pairs_symmetry, 0);
      element_encode[ct3] = 1.0;
      encode_dict[it1->GetString() + '-' + it2->GetString()] = element_encode;
      ++ct3;
    }
  }

  return encode_dict;
}

// Returns forward and backward sorted lattice lists
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMMM(
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
    if (lattice.GetId() == lattice_id_jump_pair.first
        || lattice.GetId() == lattice_id_jump_pair.second)
      continue;
    lattice_list.push_back(lattice);
  }
  RotateLatticeVector(lattice_list,
                      GetLatticePairRotationMatrix(config, lattice_id_jump_pair));
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              return PositionCompareMMM(lhs, rhs);
            });
  return lattice_list;
}
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
    if (lattice.GetId() == lattice_id_jump_pair.first
        || lattice.GetId() == lattice_id_jump_pair.second)
      continue;
    lattice_list.push_back(lattice);
  }
  RotateLatticeVector(lattice_list,
                      GetLatticePairRotationMatrix(config, lattice_id_jump_pair));
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
  /// pairs
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
      switch (FindDistanceLabelBetweenLattice(id1, id2, config)) {
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
  /// singlets
  std::vector<Singlet_MMM_t> singlet_vector;
  /// pairs
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
      switch (FindDistanceLabelBetweenLattice(id1, id2, config)) {
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
} // namespace pred
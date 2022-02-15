#include "ClusterExpansionPredictor.h"

namespace pred {
using Singlet_Periodic_t = LatticeClusterPeriodic<1>;
using Pair_Periodic_t = LatticeClusterPeriodic<2>;
using Triplet_Periodic_t = LatticeClusterPeriodic<3>;
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
template<size_t DataSize>
static bool IsClusterSmallerSymmetrically(const LatticeClusterPeriodic<DataSize> &lhs,
                                          const LatticeClusterPeriodic<DataSize> &rhs) {
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
std::vector<cfg::Lattice> GetSortedLatticeVectorPeriodic(
    const cfg::Config &config, size_t lattice_id) {
  constexpr size_t kNumOfSites = 42;
  auto lattice_id_hashset =
      GetNeighborsLatticeIdSetOfLattice(config, lattice_id);
  const auto move_distance = Vector_t{0.5, 0.5, 0.5} -
      config.GetLatticeVector()[lattice_id].GetRelativePosition();
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
              return LatticeSortCompare(lhs, rhs);
            });
  return lattice_list;
}

template<size_t DataSize>
static void GetAverageParametersMappingFromLatticeClusterVectorHelper(
    std::vector<LatticeClusterPeriodic<DataSize> > &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  std::vector<std::vector<size_t> > cluster_index_vector;
  for (const auto &cluster: cluster_vector) {
    auto cluster_index = cluster.GetIndexVector();
    cluster_index_vector.push_back(cluster_index);
  }
  cluster_mapping.push_back(cluster_index_vector);
}
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingPeriodic(
    const cfg::Config &config) {
  const size_t lattice_id = 0;
  const auto lattice_vector = GetSortedLatticeVectorPeriodic(config, lattice_id);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  /// singlets
  std::vector<Singlet_Periodic_t> singlet_vector;
  for (auto it1 = lattice_vector.cbegin(); it1 < lattice_vector.cend(); ++it1) {
    const auto index1 = static_cast<size_t>(std::distance(lattice_vector.begin(), it1));
    cfg::Lattice lattice1(*it1);
    lattice1.SetId(index1);
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::move(singlet_vector), cluster_mapping);
  /// pairs
  std::unordered_set<Pair_Periodic_t, boost::hash<Pair_Periodic_t> > first_pair_set;
  std::unordered_set<Pair_Periodic_t, boost::hash<Pair_Periodic_t> > second_pair_set;
  std::unordered_set<Pair_Periodic_t, boost::hash<Pair_Periodic_t> > third_pair_set;
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
      std::vector<Pair_Periodic_t>(first_pair_set.begin(),
                                   first_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Pair_Periodic_t>(second_pair_set.begin(),
                                   second_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Pair_Periodic_t>(third_pair_set.begin(),
                                   third_pair_set.end()), cluster_mapping);
  /// triplets
  std::unordered_set<Triplet_Periodic_t, boost::hash<Triplet_Periodic_t> >
      first_first_first_triplets_set;
  std::unordered_set<Triplet_Periodic_t, boost::hash<Triplet_Periodic_t> >
      first_first_second_triplets_set;
  std::unordered_set<Triplet_Periodic_t, boost::hash<Triplet_Periodic_t> >
      first_first_third_triplets_set;
  std::unordered_set<Triplet_Periodic_t, boost::hash<Triplet_Periodic_t> >
      first_second_third_triplets_set;
  std::unordered_set<Triplet_Periodic_t, boost::hash<Triplet_Periodic_t> >
      first_third_third_triplets_set;
  std::unordered_set<Triplet_Periodic_t, boost::hash<Triplet_Periodic_t> >
      second_third_third_triplets_set;
  std::unordered_set<Triplet_Periodic_t, boost::hash<Triplet_Periodic_t> >
      third_third_third_triplets_set;
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
          first_first_first_triplets_set.emplace(
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
          first_first_second_triplets_set.emplace(
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
          first_first_third_triplets_set.emplace(
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
          first_second_third_triplets_set.emplace(
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
          first_third_third_triplets_set.emplace(
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
          second_third_third_triplets_set.emplace(
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
          third_third_third_triplets_set.emplace(
              std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
        }
      }
    }
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Triplet_Periodic_t>(first_first_first_triplets_set.begin(),
                                      first_first_first_triplets_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Triplet_Periodic_t>(first_first_second_triplets_set.begin(),
                                      first_first_second_triplets_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Triplet_Periodic_t>(first_first_third_triplets_set.begin(),
                                      first_first_third_triplets_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Triplet_Periodic_t>(first_second_third_triplets_set.begin(),
                                      first_second_third_triplets_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Triplet_Periodic_t>(first_third_third_triplets_set.begin(),
                                      first_third_third_triplets_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Triplet_Periodic_t>(second_third_third_triplets_set.begin(),
                                      second_third_third_triplets_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelper(
      std::vector<Triplet_Periodic_t>(third_third_third_triplets_set.begin(),
                                      third_third_third_triplets_set.end()), cluster_mapping);
  return cluster_mapping;
}
std::vector<double> GetOneHotParametersFromMapPeriodic(
    const cfg::Config &config,
    const std::set<Element> &type_set,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  auto one_hot_encode_hashmap = GetOneHotEncodeHashmap(type_set);

  auto sorted_lattice_vector_periodic = GetSortedLatticeVectorPeriodic(
      config, cfg::GetVacancyLatticeIndex(config));
  // for (auto o: sorted_lattice_vector_periodic) {
  //   std::cerr << config.GetElementAtLatticeId(o.GetId()).GetString() << ' ';
  // }
  const size_t num_of_elements = type_set.size();
  std::vector<double> res_encode;
  res_encode.reserve(354);
  for (const auto &cluster_vector: cluster_mapping) {
    std::vector<double> sum_of_list(
        static_cast<size_t>(std::pow(num_of_elements, cluster_vector[0].size())), 0);
    for (const auto &cluster: cluster_vector) {
      std::string cluster_type;
      for (auto index: cluster) {
        cluster_type += config.GetElementAtLatticeId(
            sorted_lattice_vector_periodic[index].GetId()).GetString();
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
} // namespace pred
/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/3/23 11:31 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/30/23 3:13 PM                                                          *
 **************************************************************************************************/

/*! \file  ClusterExpansion.cpp
 *  \brief File for the implementation of ClusterExpansion methods.
 */

#include "ClusterExpansion.h"

/*! \brief Add neighboring sites to the existing clusters list
 *  \param reference_config : The configuration the code works on
 *  \param max_bond_order   : The cutoff bond distance order applied
 *  \param old_clusters     : The input cluster list
 *  \return                 : The output cluster list
 */
static std::vector<std::vector<size_t>> AddOneSiteToExistingClusterHelper(
    const Config &reference_config,
    const size_t max_bond_order,
    const std::vector<std::vector<size_t>> &old_clusters) {

  // Get the neighbor lists with at least m neighbors for each atom
  const auto &neighbor_lists = reference_config.GetNeighborLists();
  std::vector<std::vector<size_t>> new_clusters;
// #pragma omp parallel default(none) shared(new_clusters, old_clusters, neighbor_lists, max_bond_order, reference_config)
//   {
// #pragma omp for schedule(dynamic)
  for (const auto &old_cluster : old_clusters) {
    std::set<size_t> neighbors{};
    for (size_t m = 0; m < max_bond_order; m++) {
      // retrieve m-th nearest neighbors of i-th atom
      for (auto lattice_id : old_cluster) {
        neighbors.insert(neighbor_lists[m][lattice_id].begin(), neighbor_lists[m][lattice_id].end());
      }
    }
    // remove the sites already in the cluster
    for (auto lattice_id : old_cluster) {
      neighbors.erase(lattice_id);
    }
    // add new sites to the cluster
    for (auto new_lattice_id : neighbors) {
      std::vector<size_t> new_cluster{old_cluster};
      if (std::all_of(old_cluster.begin(), old_cluster.end(), [&](size_t old_lattice_id) {
        auto distance_order = reference_config.GetDistanceOrder(old_lattice_id, new_lattice_id);
        return old_lattice_id < new_lattice_id && distance_order <= max_bond_order /*&& distance_order > 0*/;
      })) {
        new_cluster.push_back(new_lattice_id);
// #pragma omp critical
//           {
        new_clusters.push_back(new_cluster);
//           }
//         }
      }
    }
  }
  return new_clusters;
}

LatticeClusterType IndentifyLatticeClusterType(const Config &reference_config,
                                               const std::vector<size_t> &cluster) {
  std::vector<size_t> bond_order_vector{};
  for (auto it1 = cluster.begin(); it1 != cluster.end(); it1++) {
    for (auto it2 = it1 + 1; it2 != cluster.end(); it2++) {
      bond_order_vector.push_back(reference_config.GetDistanceOrder(*it1, *it2));
    }
  }
  return LatticeClusterType{cluster.size(), bond_order_vector};
}

AtomClusterType IndentifyAtomClusterType(const Config &reference_config,
                                         const std::vector<size_t> &cluster) {
  std::vector<Element> element_vector{};
  for (const auto lattice_id : cluster) {
    element_vector.push_back(reference_config.GetElementOfLattice(lattice_id));
  }
  return AtomClusterType{element_vector};
}

ClusterType IndentifyClusterType(const Config &reference_config,
                                 const std::vector<size_t> &cluster) {
  return ClusterType{IndentifyAtomClusterType(reference_config, cluster),
                     IndentifyLatticeClusterType(reference_config, cluster)};
}

std::set<LatticeClusterType> InitializeLatticeClusterTypeSet(const Config &reference_config,
                                                             size_t max_cluster_size,
                                                             size_t max_bond_order) {
  std::set<LatticeClusterType> initialized_lattice_cluster_set;
  // Loop over lattice_type. Use one point to find all connected sites.
  std::vector<std::vector<size_t>> cluster_list{{}, {0}};
  for (size_t i = 0; i < max_cluster_size; i++) {
    if (i > 0) {
      cluster_list = AddOneSiteToExistingClusterHelper(reference_config, max_bond_order, cluster_list);
    }
    for (auto &cluster : cluster_list) {
      initialized_lattice_cluster_set.insert(IndentifyLatticeClusterType(reference_config, cluster));
    }
  }
  return initialized_lattice_cluster_set;
}

std::set<AtomClusterType> InitializeAtomClusterTypeSet(const std::set<Element> &element_set, size_t max_cluster_size) {
  std::set<AtomClusterType> initialized_atom_cluster_set;
  size_t num_elements = element_set.size();
  const std::vector<Element> element_vector(element_set.begin(), element_set.end());
  for (size_t len = 0; len <= max_cluster_size; ++len) {
    std::vector<size_t> indices(len, 0);
    while (true) {
      // Add the current combination
      std::vector<Element> combination;
      for (size_t i = 0; i < len; ++i) {
        combination.emplace_back(element_vector[indices[i]]);
      }
      auto count_vacancy = std::count(combination.begin(), combination.end(), Element("X"));
      // Remove the combination if it contains more than one vacancy
      if (count_vacancy <= 1) { initialized_atom_cluster_set.emplace(combination); }
      // Move to the next combination
      int pos = static_cast<int>(len - 1);
      while (pos >= 0 && ++indices[static_cast<size_t>(pos)] == num_elements) {
        indices[static_cast<size_t>(pos)] = 0;
        --pos;
      }
      if (pos < 0) {
        break; // All combinations generated for this length
      }
    }
  }
  return initialized_atom_cluster_set;
}

std::set<ClusterType> InitializeClusterTypeSet(const Config &reference_config,
                                               const std::set<Element> &element_set,
                                               size_t max_cluster_size,
                                               size_t max_bond_order) {
  std::set<ClusterType> initialized_cluster_set;
  auto lattice_cluster_type_set = InitializeLatticeClusterTypeSet(reference_config, max_cluster_size, max_bond_order);
  auto atom_cluster_type_set = InitializeAtomClusterTypeSet(element_set, max_cluster_size);

  for (const LatticeClusterType &lattice_cluster_type : lattice_cluster_type_set) {
    for (const AtomClusterType &atom_cluster_type : atom_cluster_type_set) {
      if (lattice_cluster_type.GetSize() != atom_cluster_type.GetSize()) { continue; }
      initialized_cluster_set.emplace(atom_cluster_type, lattice_cluster_type);
    }
  }
  return initialized_cluster_set;
}

std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> FindAllLatticeClusters(
    const Config &reference_config,
    size_t max_cluster_size,
    size_t max_bond_order,
    const std::vector<size_t> &lattice_id_vector) {
  std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> lattice_cluster_hashset;
  std::vector<std::vector<size_t>> cluster_list{{}};
  if (lattice_id_vector.empty()) {
    for (size_t i = 0; i < reference_config.GetNumLattices(); ++i) {
      cluster_list.push_back({i});
    }
  } else {
    for (auto lattice_id : lattice_id_vector) {
      cluster_list.push_back({lattice_id});
    }
  }

  for (size_t i = 0; i < max_cluster_size; i++) {
    if (i > 0) {
      cluster_list = AddOneSiteToExistingClusterHelper(reference_config, max_bond_order, cluster_list);
    }
// #pragma omp parallel default(none) shared(cluster_list, reference_config, lattice_cluster_hashset)
// {
// #pragma omp for schedule(dynamic)
    for (auto &cluster : cluster_list) {
      auto type = IndentifyLatticeClusterType(reference_config, cluster);
// #pragma omp critical
//         {
      lattice_cluster_hashset.emplace(type, cluster);
// }
    }
// }
  }
  return lattice_cluster_hashset;
}

std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> CountLatticeSite(
    const Config &reference_config, const size_t max_cluster_size, const size_t max_bond_order) {
  std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> initialized_cluster_hashmap;

  // Loop over lattice_type
  std::vector<std::vector<size_t>> cluster_list;
  for (size_t i = 0; i < reference_config.GetNumLattices(); ++i) {
    cluster_list.push_back({i});
  }
  for (size_t i = 0; i < max_cluster_size; i++) {
    if (i > 0) {
      cluster_list = AddOneSiteToExistingClusterHelper(reference_config, max_bond_order, cluster_list);
    }
    for (auto &cluster : cluster_list) {
      initialized_cluster_hashmap[IndentifyLatticeClusterType(reference_config, cluster)]++;
    }
  }
  return initialized_cluster_hashmap;
}

std::unordered_map<LatticeClusterType,
                   std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>,
                   boost::hash<LatticeClusterType>> LatticeSiteHashMap(
    const Config &reference_config, const size_t max_cluster_size, const size_t max_bond_order) {
  std::unordered_map<LatticeClusterType,
                     std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>,
                     boost::hash<LatticeClusterType>>
      initialized_cluster_hashmap;

  // Loop over lattice_type
  std::vector<std::vector<size_t>> tmp;
  for (size_t i = 0; i < reference_config.GetNumLattices(); ++i) {
    tmp.push_back({i});
  }
  for (const auto &cluster : tmp) {
    auto type = IndentifyLatticeClusterType(reference_config, cluster);
    if (initialized_cluster_hashmap.find(type) == initialized_cluster_hashmap.end())
      initialized_cluster_hashmap[type] = std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>{};
    initialized_cluster_hashmap.at(type).emplace(type, cluster);
  }
  for (size_t i = 1; i < max_cluster_size; i++) {
    tmp = AddOneSiteToExistingClusterHelper(reference_config, max_bond_order, tmp);
    for (const auto &cluster : tmp) {
      auto type = IndentifyLatticeClusterType(reference_config, cluster);
      if (initialized_cluster_hashmap.find(type) == initialized_cluster_hashmap.end()) {
        initialized_cluster_hashmap[type] = std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>{};
      }
      initialized_cluster_hashmap.at(type).emplace(type, cluster);
    }
  }
  return initialized_cluster_hashmap;
}

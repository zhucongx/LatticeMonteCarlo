/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/3/23 11:31 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/27/23 10:54 PM                                                          *
 **************************************************************************************************/

/*! \file  ClusterExpansion.cpp
 *  \brief File for the ClusterExpansion class implementation.
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
        new_clusters.push_back(new_cluster);
      }
    }
  }
  return new_clusters;
}

/*! \brief Indentify the type of a lattice cluster
 *  \param reference_config : The configuration the code works on
 *  \param cluster          : The input cluster
 *  \return                 : The lattice type of the cluster
 */
static LatticeClusterType IndentifyLatticeClusterType(const Config &reference_config,
                                                      const std::vector<size_t> &cluster) {
  std::vector<size_t> bond_order_vector{};
  for (auto it1 = cluster.begin(); it1 != cluster.end(); it1++) {
    for (auto it2 = it1 + 1; it2 != cluster.end(); it2++) {
      bond_order_vector.push_back(reference_config.GetDistanceOrder(*it1, *it2));
    }
  }
  return {cluster.size(), bond_order_vector};
}

std::set<LatticeClusterType> InitializeLatticeClusterTypeSet(const Config &reference_config,
                                                             size_t max_cluster_size,
                                                             size_t max_bond_order) {

  std::set<LatticeClusterType> initialized_lattice_cluster_set;
  // Loop over lattice_type. Use one point to find all connected sites.
  std::vector<std::vector<size_t>> cluster_list{{0}};
  for (size_t i = 0; i < max_cluster_size; i++) {
    if (i > 0) {
      cluster_list = AddOneSiteToExistingClusterHelper(reference_config, max_bond_order, cluster_list);
    }
    for (auto &cluster : cluster_list) {
      initialized_lattice_cluster_set.insert(IndentifyLatticeClusterType(reference_config, cluster));
    }
  }
  initialized_lattice_cluster_set.emplace(0, std::vector<size_t>{});
  return initialized_lattice_cluster_set;
}

std::set<AtomClusterType> InitializeAtomClusterTypeSet(const std::set<Element> &element_set) {
  std::set<AtomClusterType> initialized_atom_cluster_set;
  std::vector<Element> element_vector(element_set.begin(), element_set.end());
  size_t n = element_vector.size();
  size_t num_combinations = 1 << n; // 2^n

  for (size_t i = 0; i < num_combinations; i++) {
    std::vector<Element> combination;
    for (size_t j = 0; j < n; j++) {
      if (i & (1 << j)) {
        combination.push_back(element_vector[j]);
      }
    }
    initialized_atom_cluster_set.emplace(combination);
  }
  return initialized_atom_cluster_set;
}

std::set<ClusterType> InitializeClusterTypeSet(const Config &reference_config,
                                               size_t max_cluster_size,
                                               size_t max_bond_order,
                                               const std::set<Element> &element_set) {
  std::set<ClusterType> initialized_cluster_set;
  auto lattice_cluster_type_set = InitializeLatticeClusterTypeSet(reference_config, max_cluster_size, max_bond_order);
  auto atom_cluster_type_set = InitializeAtomClusterTypeSet(element_set);

  for (const LatticeClusterType &lattice_cluster_type : lattice_cluster_type_set) {
    for (const AtomClusterType &atom_cluster_type : atom_cluster_type_set) {
      if (lattice_cluster_type.GetSize() != atom_cluster_type.GetSize()) { continue; }
      initialized_cluster_set.emplace(atom_cluster_type, lattice_cluster_type);
    }
  }
  return initialized_cluster_set;
}

std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> CountLatticeSite(
    const Config &reference_config, const size_t max_cluster_size, const size_t max_bond_order) {
  std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> initialized_cluster_hashmap;

  // Loop over lattice_type
  std::vector<std::vector<size_t>> cluster_list;
  for (size_t i = 0; i < reference_config.GetNumSites(); ++i) {
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
  for (size_t i = 0; i < reference_config.GetNumSites(); ++i) {
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
      if (initialized_cluster_hashmap.find(type) == initialized_cluster_hashmap.end())
        initialized_cluster_hashmap[type] = std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>{};
      initialized_cluster_hashmap.at(type).emplace(type, cluster);
    }
  }
  return initialized_cluster_hashmap;
}


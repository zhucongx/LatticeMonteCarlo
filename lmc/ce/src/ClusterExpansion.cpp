/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/3/23 11:31 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/5/23 12:23 AM                                                           *
 **************************************************************************************************/

/*! \file  ClusterExpansion.cpp
 *  \brief File for the ClusterExpansion class implementation.
 */

#include "ClusterExpansion.h"

static std::vector<std::vector<size_t>> AddOneSiteToExistingClusterHelper(
    const Config &reference_config,
    const size_t max_bond_order,
    const std::vector<std::vector<size_t>> &valid_clusters) {

  // Get the neighbor lists with at least m neighbors for each atom
  const auto &neighbor_lists = reference_config.GetNeighborLists();
  std::vector<std::vector<size_t>> new_clusters;
  for (const auto &cluster : valid_clusters) {
    std::set<size_t> neighbors{};
    for (size_t m = 0; m < max_bond_order; m++) {
      // retrieve m-th nearest neighbors of i-th atom
      for (auto lattice_id : cluster) {
        neighbors.insert(neighbor_lists[m][lattice_id].begin(), neighbor_lists[m][lattice_id].end());
      }
    }
    // remove the sites already in the cluster
    for (auto it = neighbors.begin(); it != neighbors.end();) {
      if (std::find(cluster.begin(), cluster.end(), *it) != cluster.end()) {
        it = neighbors.erase(it);
      } else {
        ++it;
      }
    }
    // add the new sites to the cluster
    for (auto new_lattice_id : neighbors) {
      std::vector<size_t> temp{cluster};
      if (std::all_of(cluster.begin(), cluster.end(), [&](size_t old_lattice_id) {
        return reference_config.GetDistanceOrder(old_lattice_id, new_lattice_id) <= max_bond_order
            && old_lattice_id < new_lattice_id;
      })) {
        temp.push_back(new_lattice_id);
        new_clusters.push_back(temp);
      }
    }
  }
  return new_clusters;
}
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
std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> InitializeLatticeSiteHashMap(
    const Config &reference_config, const size_t max_cluster_size, const size_t max_bond_order) {
  std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> initialized_cluster_hashmap;

  // Loop over lattice_type
  std::vector<std::vector<size_t>> tmp;
  for (size_t i = 0; i < reference_config.GetNumSites(); ++i) {
    tmp.push_back({i});
  }
  for (auto &cluster : tmp) {
    initialized_cluster_hashmap[IndentifyLatticeClusterType(reference_config, cluster)]++;
  }
  for (size_t i = 1; i < max_cluster_size; i++) {
    tmp = AddOneSiteToExistingClusterHelper(reference_config, max_bond_order, tmp);
    for (auto &cluster : tmp) {
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


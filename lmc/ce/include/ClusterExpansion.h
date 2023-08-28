/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/3/23 11:31 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/27/23 10:10 PM                                                          *
 **************************************************************************************************/

/*! \file  ClusterExpansion.h
 *  \brief File for the ClusterExpansion class definition.
 */

#ifndef LMC_CE_INCLUDE_CLUSTEREXPANSION_H_
#define LMC_CE_INCLUDE_CLUSTEREXPANSION_H_

#include <set>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>
#include "LatticeCluster.hpp"
#include "AtomCluster.hpp"
#include "ClusterType.hpp"
/*! \brief Class for defining cluster expansion Hamiltonian.
 */
class ClusterExpansion {

  std::set<Element> element_set_{};
  Eigen::VectorXd base_theta_{};
};

std::set<LatticeClusterType> InitializeLatticeClusterTypeSet(const Config &reference_config,
                                                             size_t max_cluster_size,
                                                             size_t max_bond_order);

std::set<AtomClusterType> InitializeAtomClusterTypeSet(const std::set<Element> &element_set);

std::set<ClusterType> InitializeClusterTypeSet(const Config &reference_config,
                                               size_t max_cluster_size,
                                               size_t max_bond_order,
                                               const std::set<Element> &element_set);

std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> CountLatticeSite(
    const Config &reference_config, size_t max_cluster_size, size_t max_bond_order);
std::unordered_map<LatticeClusterType,
                   std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>,
                   boost::hash<LatticeClusterType>> LatticeSiteHashMap(
    const Config &reference_config, size_t max_cluster_size, size_t max_bond_order);

#endif //LMC_CE_INCLUDE_CLUSTEREXPANSION_H_

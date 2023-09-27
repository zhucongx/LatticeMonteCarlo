/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/3/23 11:31 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/27/23 2:15 PM                                                           *
 **************************************************************************************************/

/*! \file  ClusterExpansion.h
 *  \brief File for functions of ClusterExpansion methods.
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

/*! \brief Indentify the type of a lattice cluster
 *  \param reference_config : The configuration the code works on
 *  \param cluster          : The input cluster, of lattice ids
 *  \return                 : The lattice type of the cluster
 */
LatticeClusterType IndentifyLatticeClusterType(const Config &reference_config,
                                               const std::vector<size_t> &cluster);

/*! \brief Indentify the type of a atom cluster
 *  \param reference_config : The configuration the code works on
 *  \param cluster          : The input cluster, of lattice ids
 *  \return                 : The atom type of the cluster
 */
AtomClusterType IndentifyAtomClusterType(const Config &reference_config,
                                         const std::vector<size_t> &cluster);

/*! \brief Indentify the type of a cluster
 *  \param reference_config : The configuration the code works on
 *  \param cluster          : The input cluster
 *  \return                 : The type of the cluster
 */
ClusterType IndentifyClusterType(const Config &reference_config,
                                 const std::vector<size_t> &cluster);

/*! \brief Create a set that contains all available lattice cluster types
 *  \param reference_config : The configuration the code works on
 *  \param max_cluster_size : The maximum cluster size
 *  \param max_bond_order   : The cutoff bond distance order applied
 *  \return                 : A set of lattice cluster types
 */
std::set<LatticeClusterType> InitializeLatticeClusterTypeSet(const Config &reference_config,
                                                             size_t max_cluster_size,
                                                             size_t max_bond_order);
/*! \brief Create a set that contains all available atom cluster types
 *  \param element_set : The set of elements
 *  \return            : A set of atom cluster types
 */
std::set<AtomClusterType> InitializeAtomClusterTypeSet(const std::set<Element> &element_set, size_t max_cluster_size);

/*! \brief Create a set that contains all available cluster types
 *  \param reference_config : The configuration the code works on
 *  \param element_set      : The set of elements
 *  \param max_cluster_size : The maximum cluster size
 *  \param max_bond_order   : The cutoff bond distance order applied
 *  \return                 : A set of cluster types
 */
std::set<ClusterType> InitializeClusterTypeSet(const Config &reference_config,
                                               const std::set<Element> &element_set,
                                               size_t max_cluster_size,
                                               size_t max_bond_order);

/* \brief Create a hashset that contains all available lattice clusters related to the given lattice id vector
 * \param reference_config : The configuration the code works on
 * \param max_cluster_size : The maximum cluster size
 * \param max_bond_order   : The cutoff bond distance order applied
 * \param lattice_id_vector: The vector of lattice ids
 * \return                 : A hashset of lattice clusters
 */
std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> FindAllLatticeClusters(
    const Config &reference_config,
    size_t max_cluster_size,
    size_t max_bond_order,
    const std::vector<size_t> &lattice_id_vector);

/*! \brief Count how many lattice clusters of each type appears in the configuration
 *  \param reference_config : The configuration the code works on
 *  \param max_cluster_size : The maximum cluster size
 *  \param max_bond_order   : The cutoff bond distance order applied
 *  \return                 : A hashmap that maps the lattice cluster type to its count
 */
std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> CountLatticeSite(
    const Config &reference_config, size_t max_cluster_size, size_t max_bond_order);

/*! \brief Map the lattice sites to their corresponding lattice cluster types
 *  \param reference_config : The configuration the code works on
 *  \param max_cluster_size : The maximum cluster size
 *  \param max_bond_order   : The cutoff bond distance order applied
 *  \return                 : A hashmap that maps the lattice cluster type to its clusters set
 */
std::unordered_map<LatticeClusterType,
                   std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>,
                   boost::hash<LatticeClusterType>> LatticeSiteHashMap(
    const Config &reference_config, size_t max_cluster_size, size_t max_bond_order);

#endif //LMC_CE_INCLUDE_CLUSTEREXPANSION_H_

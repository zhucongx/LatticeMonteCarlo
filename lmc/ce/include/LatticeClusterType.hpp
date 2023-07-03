/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/2/23 3:00 PM                                                                          *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/2/23 11:06 PM                                                           *
 **************************************************************************************************/

/*! \file  LatticeClusterType.h
 *  \brief File for the LatticeClusterType class definition.
 */

#ifndef LMC_CE_INCLUDE_LATTICECLUSTERTYPE_HPP_
#define LMC_CE_INCLUDE_LATTICECLUSTERTYPE_HPP_

#include <string>
#include <vector>
#include <type_traits>
#include <boost/functional/hash.hpp>

/*! \brief Class for defining type of lattice clusters for cluster expansion.
 */
class LatticeClusterType {
 public:
  /*! \brief Default constructor for LatticeClusterType.
   */
  LatticeClusterType() = default;

  /*! \brief Constructor for setting up the LatticeClusterType.
   *  \param size                 : The size of the cluster.
   *  \param bond_distance_vector : The bond distance vector of the cluster.
   */
  LatticeClusterType(size_t size, std::vector<size_t> bond_distance_vector)
      : size_(size), bond_distance_vector_(std::move(bond_distance_vector)) {
    if (size_ * (size_ - 1) != 2 * bond_distance_vector_.size()) {
      throw std::invalid_argument("The size of the cluster is not compatible with the bond distance vector.");
    }
    std::sort(bond_distance_vector_.begin(), bond_distance_vector_.end());
  }

  /*! \brief Constructor for setting up LatticeClusterType.
   *  \param size                 : The size of the cluster.
   *  \param bond_distance_index  : The bond distance index.
   */
  template<typename ... Id>
  explicit LatticeClusterType(size_t size, Id &&... bond_distance_index)
      : size_(size), bond_distance_vector_{std::forward<Id>(bond_distance_index)...} {
    if (size_ * (size_ - 1) != 2 * bond_distance_vector_.size()) {
      throw std::invalid_argument("The size of the cluster is not compatible with the bond distance vector.");
    }
    std::sort(bond_distance_vector_.begin(), bond_distance_vector_.end());
  }

  /*! \brief 'equals' operator.
   *  \param lhs : The left hand side ClusterType.
   *  \param rhs : The right hand side ClusterType.
   *  \return    : lhs == rhs
   */
  friend bool operator==(const LatticeClusterType &lhs, const LatticeClusterType &rhs) {
    if (lhs.bond_distance_vector_.size() != rhs.bond_distance_vector_.size()) { return false; }
    for (size_t i = 0; i < lhs.bond_distance_vector_.size(); ++i) {
      if (lhs.bond_distance_vector_[i] != rhs.bond_distance_vector_[i]) { return false; }
    }
    return true;
  }

  /*! \brief Generate the hash of a ClusterType.
   *  \param cluster_type : The ClusterType to be hashed.
   *  \returns            : The hash value of the ClusterType.
   */
  friend std::size_t hash_value(const LatticeClusterType &cluster_type) {
    std::size_t seed = 0;
    boost::hash_combine(seed, cluster_type.size_);
    for (const auto &bond_distance : cluster_type.bond_distance_vector_) {
      boost::hash_combine(seed, bond_distance);
    }
    return seed;
  }

 private:
  /// The size of the cluster.
  size_t size_{};
  /// The bond distance vector of the cluster. The digits, i, in the vector represent the i^th bond distance
  std::vector<size_t> bond_distance_vector_{};
};

#endif //LMC_CE_INCLUDE_LATTICECLUSTERTYPE_HPP_

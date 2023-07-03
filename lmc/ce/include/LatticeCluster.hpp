/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 12/6/21 8:55 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/2/23 10:43 PM                                                           *
 **************************************************************************************************/

/*! \file  LatticeCluster.h
 *  \brief File for the LatticeCluster class definition.
 */

#ifndef LMC_CE_INCLUDE_LATTICECLUSTER_HPP_
#define LMC_CE_INCLUDE_LATTICECLUSTER_HPP_

#include <utility>
#include <vector>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include "LatticeClusterType.hpp"

/*! \brief Class for defining a cluster of lattice sites.
 */
class LatticeCluster {
 public:
  /*! \brief Default constructor for LatticeCluster.
   */
  LatticeCluster() = default;

  /*! \brief Constructor for setting up the cluster of lattice sites.
   *  \param lattice_cluster_type : The cluster type.
   *  \param lattice_id_vector    : The lattice id vector of the cluster
   */
  LatticeCluster(LatticeClusterType lattice_cluster_type, std::vector<size_t> lattice_id_vector)
      : cluster_type_(std::move(lattice_cluster_type)), lattice_id_vector_(std::move(lattice_id_vector)) {
    std::sort(lattice_id_vector_.begin(), lattice_id_vector_.end());
  }

  /*! \brief Constructor for setting up the cluster type.
   *  \param cluster_type : The cluster type.
   *  \param lattice_id   : The bond distance index.
   */
  template<typename ... Id>
  explicit LatticeCluster(LatticeClusterType cluster_type, Id &&... lattice_id)
      : cluster_type_(std::move(cluster_type)), lattice_id_vector_{std::forward<Id>(lattice_id)...} {
    std::sort(lattice_id_vector_.begin(), lattice_id_vector_.end());
  }

  /*! \brief Default destructor for LatticeCluster.
   */
  virtual ~LatticeCluster() = default;

  /*! \brief Query for the vector of id index of lattice sites.
   *  \return : The vector of id index of lattice sites.
   */
  [[nodiscard]] const std::vector<size_t> &GetLatticeIdVector() const {
    return lattice_id_vector_;
  }

  /*! \brief 'equals' operator.
   *  \param lhs : The left hand side LatticeCluster.
   *  \param rhs : The right hand side LatticeCluster.
   *  \return    : lhs == rhs
   */
  friend bool operator==(const LatticeCluster &lhs, const LatticeCluster &rhs) {
    // if (lhs.cluster_type_ != rhs.cluster_type_) { return false; } // not necessary
    if (lhs.lattice_id_vector_.size() != rhs.lattice_id_vector_.size()) { return false; }
    for (size_t i = 0; i < lhs.lattice_id_vector_.size(); ++i) {
      if (lhs.lattice_id_vector_[i] != rhs.lattice_id_vector_[i]) { return false; }
    }
    return true;
  }

  /*! \brief Generate the hash of a LatticeCluster.
   *  \param cluster : The LatticeCluster to be hashed.
   *  \returns       : The hash value of the LatticeCluster.
   */
  friend size_t hash_value(const LatticeCluster &cluster) {
    size_t seed = 0;
    for (size_t i : cluster.lattice_id_vector_) {
      boost::hash_combine(seed, i);
    }
    return seed;
  }

 protected:
  /// The type of the cluster.
  LatticeClusterType cluster_type_;
  /// The vector of id index of lattice sites.
  std::vector<size_t> lattice_id_vector_;

};

#endif //LMC_CE_INCLUDE_LATTICECLUSTER_HPP_

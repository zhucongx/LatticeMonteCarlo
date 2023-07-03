/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 2/14/22 8:45 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/2/23 10:54 PM                                                           *
 **************************************************************************************************/

/*! \file  AtomCluster.h
 *  \brief File for the AtomCluster class definition.
 */

#ifndef LMC_CE_INCLUDE_ATOMCLUSTER_HPP_
#define LMC_CE_INCLUDE_ATOMCLUSTER_HPP_

#include <utility>
#include <vector>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include "AtomClusterType.hpp"

/*! \brief Class for defining a cluster of atoms.
 */
class AtomCluster {
 public:
  /*! \brief Default constructor for AtomCluster.
   */
  AtomCluster() = default;

  /*! \brief Constructor for setting up the cluster of atoms.
   *  \param cluster_type   : The cluster type.
   *  \param atom_id_vector : The atom id vector of the cluster
   */
  AtomCluster(AtomClusterType cluster_type, std::vector<size_t> atom_id_vector)
      : cluster_type_(std::move(cluster_type)), atom_id_vector_(std::move(atom_id_vector)) {
    std::sort(atom_id_vector_.begin(), atom_id_vector_.end());
  }

  /*! \brief Constructor for setting up the cluster type.
   *  \param cluster_type : The cluster type.
   *  \param atom_id      : The bond distance index.
   */
  template<typename ... Id>
  explicit AtomCluster(AtomClusterType cluster_type, Id &&... atom_id)
      : cluster_type_(std::move(cluster_type)), atom_id_vector_{std::forward<Id>(atom_id)...} {
    std::sort(atom_id_vector_.begin(), atom_id_vector_.end());
  }

  /*! \brief Default destructor for AtomCluster.
   */
  virtual ~AtomCluster() = default;

  /*! \brief Query for the vector of id index of atom .
   *  \return : The vector of id index of atoms.
   */
  [[nodiscard]] const std::vector<size_t> &GetAtomIdVector() const {
    return atom_id_vector_;
  }

  /*! \brief 'equals' operator.
   *  \param lhs : The left hand side AtomCluster.
   *  \param rhs : The right hand side AtomCluster.
   *  \return    : lhs == rhs
   */
  friend bool operator==(const AtomCluster &lhs, const AtomCluster &rhs) {
    // if (lhs.cluster_type_ != rhs.cluster_type_) { return false; } // not necessary
    if (lhs.atom_id_vector_.size() != rhs.atom_id_vector_.size()) { return false; }
    for (size_t i = 0; i < lhs.atom_id_vector_.size(); ++i) {
      if (lhs.atom_id_vector_[i] != rhs.atom_id_vector_[i]) { return false; }
    }
    return true;
  }

  /*! \brief Generate the hash of a AtomCluster.
   *  \param cluster : The AtomCluster to be hashed.
   *  \returns       : The hash value of the AtomCluster.
   */
  friend size_t hash_value(const AtomCluster &cluster) {
    size_t seed = 0;
    for (size_t i : cluster.atom_id_vector_) {
      boost::hash_combine(seed, i);
    }
    return seed;
  }

 protected:
  /// The type of the cluster.
  AtomClusterType cluster_type_;
  /// The vector of id index of atoms.
  std::vector<size_t> atom_id_vector_;

};

#endif //LMC_CE_INCLUDE_ATOMCLUSTER_HPP_

/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/2/23 3:00 PM                                                                          *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/27/23 11:08 PM                                                          *
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
   *  \param size              : The size of the cluster.
   *  \param bond_order_vector : The bond distance vector of the cluster.
   */
  LatticeClusterType(size_t size, std::vector<size_t> bond_order_vector)
      : size_(size), bond_order_vector_(std::move(bond_order_vector)) {
    if (size_ * (size_ - 1) != 2 * bond_order_vector_.size()) {
      throw std::invalid_argument("The size of the cluster is not compatible with the bond distance vector.");
    }
    std::sort(bond_order_vector_.begin(), bond_order_vector_.end());
  }

  /*! \brief Constructor for setting up LatticeClusterType.
   *  \param size       : The size of the cluster.
   *  \param bond_order : The bond distance index.
   */
  template<typename ... Id>
  explicit LatticeClusterType(size_t size, Id &&... bond_order)
      : size_(size), bond_order_vector_{std::forward<Id>(bond_order)...} {
    if (size_ * (size_ - 1) != 2 * bond_order_vector_.size()) {
      throw std::invalid_argument("The size of the cluster is not compatible with the bond distance vector.");
    }
    std::sort(bond_order_vector_.begin(), bond_order_vector_.end());
  }

  /*! \brief Equality operator.
   *  \param lhs : The left hand side LatticeClusterType.
   *  \param rhs : The right hand side LatticeClusterType.
   *  \return    : lhs == rhs
   */
  friend bool operator==(const LatticeClusterType &lhs, const LatticeClusterType &rhs) {
    if (lhs.bond_order_vector_.size() != rhs.bond_order_vector_.size()) { return false; }
    for (size_t i = 0; i < lhs.bond_order_vector_.size(); ++i) {
      if (lhs.bond_order_vector_[i] != rhs.bond_order_vector_[i]) { return false; }
    }
    return true;
  }

  /*! \brief Relational operator.
   *  \param lhs : The left hand side LatticeClusterType.
   *  \param rhs : The right hand side LatticeClusterType.
   *  \return    : lhs < rhs
   */
  friend bool operator<(const LatticeClusterType &lhs, const LatticeClusterType &rhs) {
    if (lhs.size_ < rhs.size_) { return true; }
    if (rhs.size_ < lhs.size_) { return false; }
    for (size_t i = 0; i < lhs.bond_order_vector_.size(); ++i) {
      if (lhs.bond_order_vector_[i] < rhs.bond_order_vector_[i]) { return true; }
      if (rhs.bond_order_vector_[i] < lhs.bond_order_vector_[i]) { return false; }
    }
    return false;
  }

  /*! \brief Generate the hash of a LatticeClusterType.
   *  \param cluster_type : The LatticeClusterType to be hashed.
   *  \returns            : The hash value of the LatticeClusterType.
   */
  friend std::size_t hash_value(const LatticeClusterType &cluster_type) {
    std::size_t seed = 0;
    boost::hash_combine(seed, cluster_type.size_);
    for (const auto &bond_order : cluster_type.bond_order_vector_) {
      boost::hash_combine(seed, bond_order);
    }
    return seed;
  }

  /*! \brief stream operator output for LatticeClusterType.
   *  \param os           : The output stream.
   *  \param cluster_type : The LatticeClusterType to be streamed.
   *  \return             : The output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, const LatticeClusterType &cluster_type) {
    os << cluster_type.size_;
    for (auto bond_order : cluster_type.bond_order_vector_) {
      os << '-' << bond_order;
    }
    return os;
  }

  /*! \brief Query for the size of the cluster.
   *  \return : The size of the cluster.
   */
  [[nodiscard]] size_t GetSize() const {
    return size_;
  }

 private:
  /// The size of the cluster.
  size_t size_{};
  /// The bond distance vector of the cluster. The digits, i, in the vector represent the i^th bond distance
  std::vector<size_t> bond_order_vector_{};
};

#endif //LMC_CE_INCLUDE_LATTICECLUSTERTYPE_HPP_

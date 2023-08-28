/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/2/23 10:43 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/27/23 11:08 PM                                                          *
 **************************************************************************************************/

/*! \file  AtomClusterType.h
 *  \brief File for the AtomClusterType class definition.
 */

#ifndef LMC_CE_INCLUDE_ATOMCLUSTERTYPE_HPP_
#define LMC_CE_INCLUDE_ATOMCLUSTERTYPE_HPP_

#include <string>
#include <vector>
#include <type_traits>
#include <boost/functional/hash.hpp>
#include <Element.hpp>

/*! \brief Class for defining type of lattice clusters for cluster expansion.
 */
class AtomClusterType {
 public:
  /*! \brief Default constructor for AtomClusterType.
   */
  AtomClusterType() = default;

  /*! \brief Constructor for setting up the AtomClusterType.
   *  \param element_vector : The vector of the cluster.
   */
  explicit AtomClusterType(std::vector<Element> element_vector)
      : element_vector_(std::move(element_vector)) {
    std::sort(element_vector_.begin(), element_vector_.end());
  }

  /*! \brief Constructor for setting up AtomClusterType.
   *  \param size                 : The size of the cluster.
   *  \param bond_distance_index  : The bond distance index.
   */
  template<typename ... Id>
  explicit AtomClusterType(Id &&... element)
      : element_vector_{std::forward<Id>(element)...} {
    std::sort(element_vector_.begin(), element_vector_.end());
  }

  /*! \brief Equality operator.
   *  \param lhs : The left hand side AtomClusterType.
   *  \param rhs : The right hand side AtomClusterType.
   *  \return    : lhs == rhs
   */
  friend bool operator==(const AtomClusterType &lhs, const AtomClusterType &rhs) {
    if (lhs.element_vector_.size() != rhs.element_vector_.size()) { return false; }
    for (size_t i = 0; i < lhs.element_vector_.size(); ++i) {
      if (lhs.element_vector_[i] != rhs.element_vector_[i]) { return false; }
    }
    return true;
  }

  /*! \brief Relational operator.
   *  \param lhs : The left hand side AtomClusterType.
   *  \param rhs : The right hand side AtomClusterType.
   *  \return    : lhs < rhs
   */
  friend bool operator<(const AtomClusterType &lhs, const AtomClusterType &rhs) {
    if (lhs.element_vector_.size() < rhs.element_vector_.size()) { return true; }
    if (rhs.element_vector_.size() < lhs.element_vector_.size()) { return false; }
    for (size_t i = 0; i < lhs.element_vector_.size(); ++i) {
      if (lhs.element_vector_[i] < rhs.element_vector_[i]) { return true; }
      if (rhs.element_vector_[i] < lhs.element_vector_[i]) { return false; }
    }
    return false;
  }

  /*! \brief Generate the hash of a AtomClusterType.
   *  \param cluster_type : The ClusterType to be hashed.
   *  \returns            : The hash value of the AtomClusterType.
   */
  friend std::size_t hash_value(const AtomClusterType &cluster_type) {
    std::size_t seed = 0;
    for (auto element : cluster_type.element_vector_) {
      boost::hash_combine(seed, element);
    }
    return seed;
  }

  /*! \brief stream operator output for AtomClusterType.
   *  \param os           : The output stream.
   *  \param cluster_type : The AtomClusterType to be streamed.
   *  \return             : The output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, const AtomClusterType &cluster_type) {
    if (cluster_type.element_vector_.empty()) {
      os << "∅";
    }

    for (size_t i = 0; i < cluster_type.element_vector_.size(); ++i) {
      os << cluster_type.element_vector_[i];
      if (i < cluster_type.element_vector_.size() - 1) {
        os << '-';
      }
    }
    return os;
  }

  /*! \brief Query for the size of the cluster.
   *  \return : The size of the cluster.
   */
  [[nodiscard]] size_t GetSize() const {
    return element_vector_.size();
  }

 private:
  /// The chemical elements of the cluster.
  std::vector<Element> element_vector_{};
};

#endif //LMC_CE_INCLUDE_ATOMCLUSTERTYPE_HPP_

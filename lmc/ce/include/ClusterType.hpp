/**************************************************************************************************
 * Copyright (c) 2021-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 12/6/21 8:55 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/27/23 10:51 PM                                                          *
 **************************************************************************************************/

/*! \file  Cluster.h
 *  \brief File for the Cluster class definition.
 */

#ifndef LMC_CFG_INCLUDE_CLUSTERTYPE_HPP_
#define LMC_CFG_INCLUDE_CLUSTERTYPE_HPP_

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <type_traits>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "AtomClusterType.hpp"
#include "LatticeClusterType.hpp"

/*! \brief Class for defining a cluster of atoms and site positions.
 */
class ClusterType {
 public:
  ClusterType() = default;

  ClusterType(AtomClusterType atom_cluster_type, LatticeClusterType lattice_cluster_type)
      : atom_cluster_type_(std::move(atom_cluster_type)), lattice_cluster_type_(std::move(lattice_cluster_type)) {
    if (lattice_cluster_type_.GetSize() != atom_cluster_type_.GetSize()) {
      throw std::invalid_argument("The size of the lattice cluster and the atom cluster are not compatible.");
    }
  }

  friend bool operator==(const ClusterType &lhs, const ClusterType &rhs) {
    return lhs.atom_cluster_type_ == rhs.atom_cluster_type_ &&
        lhs.lattice_cluster_type_ == rhs.lattice_cluster_type_;
  }
  friend bool operator<(const ClusterType &lhs, const ClusterType &rhs) {
    if (lhs.lattice_cluster_type_ < rhs.lattice_cluster_type_)
      return true;
    if (rhs.lattice_cluster_type_ < lhs.lattice_cluster_type_)
      return false;
    return lhs.atom_cluster_type_ < rhs.atom_cluster_type_;
  }
  friend std::size_t hash_value(const ClusterType &cluster_type) {
    std::size_t seed = 0;
    boost::hash_combine(seed, cluster_type.atom_cluster_type_);
    boost::hash_combine(seed, cluster_type.lattice_cluster_type_);
    return seed;
  }

  /*! \brief stream operator output for ClusterType.
   *  \param os           : The output stream.
   *  \param cluster_type : The AtomClusterType to be streamed.
   *  \return             : The output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, const ClusterType &cluster_type) {
    os << cluster_type.lattice_cluster_type_ << " " << cluster_type.atom_cluster_type_;
    return os;
  }

 private:

  AtomClusterType atom_cluster_type_;
  LatticeClusterType lattice_cluster_type_;
};

// inline bool PositionCompareState(const cfg::Lattice &lhs,
//                                  const cfg::Lattice &rhs) {
//   const auto &relative_position_lhs = lhs.GetRelativePosition();
//   const auto &relative_position_rhs = rhs.GetRelativePosition();
//   const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
//   const double diff_y = relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
//   const double diff_z = relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
//   if (diff_x < -kEpsilon) { return true; }
//   if (diff_x > kEpsilon) { return false; }
//   if (diff_y < -kEpsilon) { return true; }
//   if (diff_y > kEpsilon) { return false; }
//   return diff_z < -kEpsilon;
// }
//
// inline bool GroupCompareMMM(const cfg::Lattice &lhs,
//                             const cfg::Lattice &rhs) {
//   const auto &relative_position_lhs = lhs.GetRelativePosition();
//   const auto &relative_position_rhs = rhs.GetRelativePosition();
//
//   const double
//       diff_norm = Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
//   if (diff_norm < -kEpsilon)
//     return true;
//   if (diff_norm > kEpsilon)
//     return false;
//   const double diff_x_sym = std::abs(relative_position_lhs[kXDimension] - 0.5)
//       - std::abs(relative_position_rhs[kXDimension] - 0.5);
//   return diff_x_sym < -kEpsilon;
// }
// inline bool GroupCompareMM2(const cfg::Lattice &lhs,
//                             const cfg::Lattice &rhs) {
//   const auto &relative_position_lhs = lhs.GetRelativePosition();
//   const auto &relative_position_rhs = rhs.GetRelativePosition();
//
//   const double
//       diff_norm = Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
//   if (diff_norm < -kEpsilon)
//     return true;
//   if (diff_norm > kEpsilon)
//     return false;
//   const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
//   return diff_x < -kEpsilon;
// }
//
// inline bool PositionCompareMMM(const cfg::Lattice &lhs,
//                                const cfg::Lattice &rhs) {
//   const auto &relative_position_lhs = lhs.GetRelativePosition();
//   const auto &relative_position_rhs = rhs.GetRelativePosition();
//
//   const double diff_norm =
//       Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
//   if (diff_norm < -kEpsilon) { return true; }
//   if (diff_norm > kEpsilon) { return false; }
//   const double diff_x_sym = std::abs(relative_position_lhs[kXDimension] - 0.5)
//       - std::abs(relative_position_rhs[kXDimension] - 0.5);
//   if (diff_x_sym < -kEpsilon) { return true; }
//   if (diff_x_sym > kEpsilon) { return false; }
//   const double diff_x =
//       relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
//   if (diff_x < -kEpsilon) { return true; }
//   if (diff_x > kEpsilon) { return false; }
//   const double diff_y =
//       relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
//   if (diff_y < -kEpsilon) { return true; }
//   if (diff_y > kEpsilon) { return false; }
//   const double diff_z =
//       relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
//   if (diff_z < -kEpsilon) { return true; }
//   if (diff_z > kEpsilon) { return false; }
//   return false;
// }
// inline bool PositionCompareMM2(const cfg::Lattice &lhs,
//                                const cfg::Lattice &rhs) {
//   const auto &relative_position_lhs = lhs.GetRelativePosition();
//   const auto &relative_position_rhs = rhs.GetRelativePosition();
//   const double diff_norm =
//       Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
//   if (diff_norm < -kEpsilon) { return true; }
//   if (diff_norm > kEpsilon) { return false; }
//   const double diff_x =
//       relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
//   if (diff_x < -kEpsilon) { return true; }
//   if (diff_x > kEpsilon) { return false; }
//   const double diff_y =
//       relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
//   if (diff_y < -kEpsilon) { return true; }
//   if (diff_y > kEpsilon) { return false; }
//   const double diff_z =
//       relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
//   if (diff_z < -kEpsilon) { return true; }
//   if (diff_z > kEpsilon) { return false; }
//   return false;
// }
// template<size_t DataSize>
// class LatticeClusterMMM : public LatticeCluster<DataSize> {
//  public:
//   explicit LatticeClusterMMM(const std::array<cfg::Lattice, DataSize> &lattice_array)
//       : LatticeCluster<DataSize>(lattice_array) {
//     Sort();
//     this->symmetry_label_ = FindSymmetryLabel();
//   }
//   ~LatticeClusterMMM() override = default;
//  private:
//   void Sort() {
//     std::sort(this->lattice_array_.begin(), this->lattice_array_.end(),
//               [](const auto &lhs, const auto &rhs) -> bool {
//                 return PositionCompareMMM(lhs, rhs);
//               });
//   }
//   bool FindSymmetryLabel() {
//     for (auto it1 = this->lattice_array_.begin(); it1 < this->lattice_array_.end(); ++it1) {
//       for (auto it2 = this->lattice_array_.begin(); it2 < it1; ++it2) {
//         if (!GroupCompareMMM(*it1, *it2) && !GroupCompareMMM(*it2, *it1)) {
//           return true;
//         }
//       }
//     }
//     return false;
//   }
// };
//
// template<size_t DataSize>
// class LatticeClusterMM2 : public LatticeCluster<DataSize> {
//  public:
//   explicit LatticeClusterMM2(const std::array<cfg::Lattice, DataSize> &lattice_array)
//       : LatticeCluster<DataSize>(lattice_array) {
//     Sort();
//     this->symmetry_label_ = FindSymmetryLabel();
//   }
//   ~LatticeClusterMM2() override = default;
//  private:
//   void Sort() {
//     std::sort(this->lattice_array_.begin(), this->lattice_array_.end(),
//               [](const auto &lhs, const auto &rhs) -> bool {
//                 return PositionCompareMM2(lhs, rhs);
//               });
//   }
//   bool FindSymmetryLabel() {
//     for (auto it1 = this->lattice_array_.begin(); it1 < this->lattice_array_.end(); ++it1) {
//       for (auto it2 = this->lattice_array_.begin(); it2 < it1; ++it2) {
//         if (!GroupCompareMM2(*it1, *it2) && !GroupCompareMM2(*it2, *it1)) {
//           return true;
//         }
//       }
//     }
//     return false;
//   }
// };
//
// template<size_t DataSize>
// inline bool IsClusterSmallerSymmetricallyMMM(const cfg::LatticeClusterMMM<DataSize> &lhs,
//                                              const cfg::LatticeClusterMMM<DataSize> &rhs) {
//   for (size_t i = 0; i < DataSize; ++i) {
//     const auto &lhs_lattice = lhs.GetLatticeAt(i);
//     const auto &rhs_lattice = rhs.GetLatticeAt(i);
//     if (GroupCompareMMM(lhs_lattice, rhs_lattice)) { return true; }
//     if (GroupCompareMMM(rhs_lattice, lhs_lattice)) { return false; }
//   }
//   // if it reaches here, it means that the clusters are same symmetrically. Returns false.
//   return false;
// }
// template<size_t DataSize>
// inline bool IsClusterSmallerSymmetricallyMM2(const cfg::LatticeClusterMM2<DataSize> &lhs,
//                                              const cfg::LatticeClusterMM2<DataSize> &rhs) {
//   for (size_t i = 0; i < DataSize; ++i) {
//     const auto &lhs_lattice = lhs.GetLatticeAt(i);
//     const auto &rhs_lattice = rhs.GetLatticeAt(i);
//     if (GroupCompareMM2(lhs_lattice, rhs_lattice)) { return true; }
//     if (GroupCompareMM2(rhs_lattice, lhs_lattice)) { return false; }
//   }
//   // if it reaches here, it means that the clusters are same symmetrically. Returns false.
//   return false;
// }

#endif //LMC_CFG_INCLUDE_CLUSTERTYPE_HPP_

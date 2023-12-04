#ifndef LMC_LMC_CFG_INCLUDE_LATTICECLUSTER_HPP_
#define LMC_LMC_CFG_INCLUDE_LATTICECLUSTER_HPP_
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <type_traits>
#include <boost/functional/hash.hpp>
#include "Lattice.hpp"
#include "Config.h"

namespace cfg {
template<size_t DataSize>
class LatticeCluster {
  public:
    /// Constructor
    explicit LatticeCluster(std::array<cfg::Lattice, DataSize> lattice_array)
        : lattice_array_(std::move(lattice_array)) {
    }
    /// Destructor
    virtual ~LatticeCluster() = default;
    /// Getter
    [[nodiscard]] const cfg::Lattice &GetLatticeAt(size_t i) const {
      return lattice_array_[i];
    }
    [[nodiscard]] bool IsSymmetryLabel() const {
      return symmetry_label_;
    }
    [[nodiscard]] std::vector<size_t> GetIndexVector() const {
      std::vector<size_t> cluster_index;
      cluster_index.reserve(DataSize);
      std::transform(lattice_array_.begin(), lattice_array_.end(),
                     std::back_inserter(cluster_index),
                     [](const auto &lattice) { return lattice.GetId(); });
      return cluster_index;
    }
    ///Operators
    friend bool operator==(const LatticeCluster<DataSize> &lhs,
                           const LatticeCluster<DataSize> &rhs) {
      for (size_t i = 0; i < DataSize; ++i) {
        if (lhs.lattice_array_[i].GetId() != rhs.lattice_array_[i].GetId())
          return false;
      }
      return true;
    }
    friend size_t hash_value(const LatticeCluster<DataSize> &cluster) {
      size_t seed = 0;
      for (size_t i = 0; i < DataSize; ++i) {
        boost::hash_combine(seed, cluster.lattice_array_[i].GetId());
      }
      return seed;
    }
  protected:
    std::array<cfg::Lattice, DataSize> lattice_array_;
    bool symmetry_label_{false};
};

inline bool PositionCompareState(const cfg::Lattice &lhs,
                                 const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();
  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  const double diff_y = relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
  const double diff_z = relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
  if (diff_x < -kEpsilon) { return true; }
  if (diff_x > kEpsilon) { return false; }
  if (diff_y < -kEpsilon) { return true; }
  if (diff_y > kEpsilon) { return false; }
  return diff_z < -kEpsilon;
}

inline bool GroupCompareMMM(const cfg::Lattice &lhs,
                            const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double
      diff_norm = Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
  if (diff_norm < -kEpsilon)
    return true;
  if (diff_norm > kEpsilon)
    return false;
  const double diff_x_sym = std::abs(relative_position_lhs[kXDimension] - 0.5)
      - std::abs(relative_position_rhs[kXDimension] - 0.5);
  return diff_x_sym < -kEpsilon;
}
inline bool GroupCompareMM2(const cfg::Lattice &lhs,
                            const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double
      diff_norm = Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
  if (diff_norm < -kEpsilon)
    return true;
  if (diff_norm > kEpsilon)
    return false;
  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  return diff_x < -kEpsilon;
}

inline bool PositionCompareMMM(const cfg::Lattice &lhs,
                               const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double diff_norm =
      Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
  if (diff_norm < -kEpsilon) { return true; }
  if (diff_norm > kEpsilon) { return false; }
  const double diff_x_sym = std::abs(relative_position_lhs[kXDimension] - 0.5)
      - std::abs(relative_position_rhs[kXDimension] - 0.5);
  if (diff_x_sym < -kEpsilon) { return true; }
  if (diff_x_sym > kEpsilon) { return false; }
  const double diff_x =
      relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  if (diff_x < -kEpsilon) { return true; }
  if (diff_x > kEpsilon) { return false; }
  const double diff_y =
      relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
  if (diff_y < -kEpsilon) { return true; }
  if (diff_y > kEpsilon) { return false; }
  const double diff_z =
      relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
  if (diff_z < -kEpsilon) { return true; }
  if (diff_z > kEpsilon) { return false; }
  return false;
}
inline bool PositionCompareMM2(const cfg::Lattice &lhs,
                               const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();
  const double diff_norm =
      Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
  if (diff_norm < -kEpsilon) { return true; }
  if (diff_norm > kEpsilon) { return false; }
  const double diff_x =
      relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  if (diff_x < -kEpsilon) { return true; }
  if (diff_x > kEpsilon) { return false; }
  const double diff_y =
      relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
  if (diff_y < -kEpsilon) { return true; }
  if (diff_y > kEpsilon) { return false; }
  const double diff_z =
      relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
  if (diff_z < -kEpsilon) { return true; }
  if (diff_z > kEpsilon) { return false; }
  return false;
}
template<size_t DataSize>
class LatticeClusterMMM : public LatticeCluster<DataSize> {
  public:
    explicit LatticeClusterMMM(const std::array<cfg::Lattice, DataSize> &lattice_array)
        : LatticeCluster<DataSize>(lattice_array) {
      Sort();
      this->symmetry_label_ = FindSymmetryLabel();
    }
    ~LatticeClusterMMM() override = default;
  private:
    void Sort() {
      std::sort(this->lattice_array_.begin(), this->lattice_array_.end(),
                [](const auto &lhs, const auto &rhs) -> bool {
                  return PositionCompareMMM(lhs, rhs);
                });
    }
    bool FindSymmetryLabel() {
      for (auto it1 = this->lattice_array_.begin(); it1 < this->lattice_array_.end(); ++it1) {
        for (auto it2 = this->lattice_array_.begin(); it2 < it1; ++it2) {
          if (!GroupCompareMMM(*it1, *it2) && !GroupCompareMMM(*it2, *it1)) {
            return true;
          }
        }
      }
      return false;
    }
};

template<size_t DataSize>
class LatticeClusterMM2 : public LatticeCluster<DataSize> {
  public:
    explicit LatticeClusterMM2(const std::array<cfg::Lattice, DataSize> &lattice_array)
        : LatticeCluster<DataSize>(lattice_array) {
      Sort();
      this->symmetry_label_ = FindSymmetryLabel();
    }
    ~LatticeClusterMM2() override = default;
  private:
    void Sort() {
      std::sort(this->lattice_array_.begin(), this->lattice_array_.end(),
                [](const auto &lhs, const auto &rhs) -> bool {
                  return PositionCompareMM2(lhs, rhs);
                });
    }
    bool FindSymmetryLabel() {
      for (auto it1 = this->lattice_array_.begin(); it1 < this->lattice_array_.end(); ++it1) {
        for (auto it2 = this->lattice_array_.begin(); it2 < it1; ++it2) {
          if (!GroupCompareMM2(*it1, *it2) && !GroupCompareMM2(*it2, *it1)) {
            return true;
          }
        }
      }
      return false;
    }
};

template<size_t DataSize>
inline bool IsClusterSmallerSymmetricallyMMM(const cfg::LatticeClusterMMM<DataSize> &lhs,
                                             const cfg::LatticeClusterMMM<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    const auto &lhs_lattice = lhs.GetLatticeAt(i);
    const auto &rhs_lattice = rhs.GetLatticeAt(i);
    if (GroupCompareMMM(lhs_lattice, rhs_lattice)) { return true; }
    if (GroupCompareMMM(rhs_lattice, lhs_lattice)) { return false; }
  }
  // if it reaches here, it means that the clusters are same symmetrically. Returns false.
  return false;
}
template<size_t DataSize>
inline bool IsClusterSmallerSymmetricallyMM2(const cfg::LatticeClusterMM2<DataSize> &lhs,
                                             const cfg::LatticeClusterMM2<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    const auto &lhs_lattice = lhs.GetLatticeAt(i);
    const auto &rhs_lattice = rhs.GetLatticeAt(i);
    if (GroupCompareMM2(lhs_lattice, rhs_lattice)) { return true; }
    if (GroupCompareMM2(rhs_lattice, lhs_lattice)) { return false; }
  }
  // if it reaches here, it means that the clusters are same symmetrically. Returns false.
  return false;
}
} // cfg

#endif //LMC_LMC_CFG_INCLUDE_LATTICECLUSTER_HPP_

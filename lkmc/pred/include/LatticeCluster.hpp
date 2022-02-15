#ifndef LKMC_LKMC_PRED_INCLUDE_LATTICECLUSTER_HPP_
#define LKMC_LKMC_PRED_INCLUDE_LATTICECLUSTER_HPP_
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <type_traits>
#include <boost/functional/hash.hpp>
#include "Lattice.hpp"
#include "Config.h"

namespace pred {
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
    virtual void Sort() = 0;
  protected:
    std::array<cfg::Lattice, DataSize> lattice_array_;
};

template<size_t DataSize>
class LatticeClusterState : public LatticeCluster<DataSize> {
  public:
    explicit LatticeClusterState(const std::array<cfg::Lattice, DataSize> &lattice_array)
        : LatticeCluster<DataSize>(lattice_array) {
      Sort();
    }
    ~LatticeClusterState() override = default;
  private:
    void Sort() override {}
};

template<size_t DataSize>
class LatticeClusterMM2 : public LatticeCluster<DataSize> {
  public:
    explicit LatticeClusterMM2(const std::array<cfg::Lattice, DataSize> &lattice_array)
        : LatticeCluster<DataSize>(lattice_array) {
      Sort();
    }
    ~LatticeClusterMM2() override = default;
  private:
    void Sort() override {
      std::sort(this->lattice_array_.begin(), this->lattice_array_.end(),
                [](const auto &lhs, const auto &rhs) {
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
                  return lhs.GetId() < rhs.GetId();
                });
    }
};

template<size_t DataSize>
class LatticeClusterMMM : public LatticeCluster<DataSize> {
  public:
    explicit LatticeClusterMMM(const std::array<cfg::Lattice, DataSize> &lattice_array)
        : LatticeCluster<DataSize>(lattice_array) {
      Sort();
    }
    ~LatticeClusterMMM() override = default;
  private:
    void Sort() override {
      std::sort(this->lattice_array_.begin(), this->lattice_array_.end(),
                [](const auto &lhs, const auto &rhs) {
                  const auto &relative_position_lhs = lhs.GetRelativePosition();
                  const auto &relative_position_rhs = rhs.GetRelativePosition();
                  const double diff_norm =
                      Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
                  if (diff_norm < -kEpsilon) { return true; }
                  if (diff_norm > kEpsilon) { return false; }
                  const double diff_x = std::abs(relative_position_lhs[kXDimension] - 0.5)
                      - std::abs(relative_position_rhs[kXDimension] - 0.5);
                  if (diff_x < -kEpsilon) { return true; }
                  if (diff_x > kEpsilon) { return false; }
                  const double diff_x2 =
                      relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
                  if (diff_x2 < -kEpsilon) { return true; }
                  if (diff_x2 > kEpsilon) { return false; }
                  const double diff_y =
                      relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
                  if (diff_y < -kEpsilon) { return true; }
                  if (diff_y > kEpsilon) { return false; }
                  const double diff_z =
                      relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
                  if (diff_z < -kEpsilon) { return true; }
                  if (diff_z > kEpsilon) { return false; }
                  return lhs.GetId() < rhs.GetId();
                });
    }
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_LATTICECLUSTER_HPP_

#ifndef LKMC_LKMC_CFG_INCLUDE_LATTICECLUSTER_HPP_
#define LKMC_LKMC_CFG_INCLUDE_LATTICECLUSTER_HPP_
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
};

// template<size_t DataSize>
// class LatticeClusterState : public LatticeCluster<DataSize> {
//   public:
//     explicit LatticeClusterState(const std::array<cfg::Lattice, DataSize> &lattice_array)
//         : LatticeCluster<DataSize>(lattice_array) {
//     }
//     ~LatticeClusterState() override = default;
//   private:
//     void Sort() override {}
// };
} // namespace cfg

#endif //LKMC_LKMC_CFG_INCLUDE_LATTICECLUSTER_HPP_

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
  public:
    std::array<cfg::Lattice, DataSize> lattice_array_;
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

template<size_t DataSize>
class LatticeClusterPeriodic : public LatticeCluster<DataSize> {
  public:
    explicit LatticeClusterPeriodic(const std::array<cfg::Lattice, DataSize> &lattice_array)
        : LatticeCluster<DataSize>(lattice_array) {
      Sort();
    }
    ~LatticeClusterPeriodic() override = default;
  private:
    void Sort() override {
      std::sort(this->lattice_array_.begin(), this->lattice_array_.end(),
                [](const auto &lhs, const auto &rhs) {
                  const auto &relative_position_lhs = lhs.GetRelativePosition();
                  const auto &relative_position_rhs = rhs.GetRelativePosition();
                  const double diff_x =
                      relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
                  const double diff_y =
                      relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
                  const double diff_z =
                      relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
                  if (diff_x < -kEpsilon) { return true; }
                  if (diff_x > kEpsilon) { return false; }
                  if (diff_y < -kEpsilon) { return true; }
                  if (diff_y > kEpsilon) { return false; }
                  if (diff_z < -kEpsilon) { return true; }
                  if (diff_z > kEpsilon) { return false; }
                  return lhs.GetId() < rhs.GetId();
                });
    }
};

inline std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<Element> &type_set) {
  size_t type_size = type_set.size();
  std::unordered_map<std::string, std::vector<double> > encode_dict;

  size_t ct1 = 0;
  for (const auto &element: type_set) {
    std::vector<double> element_encode(type_size, 0);
    element_encode[ct1] = 1.0;
    encode_dict[element.GetString()] = element_encode;
    ++ct1;
  }

  size_t num_pairs = type_size * type_size;
  size_t ct2 = 0;
  for (const auto &element1: type_set) {
    for (const auto &element2: type_set) {
      std::vector<double> element_encode(num_pairs, 0);
      element_encode[ct2] = 1.0;
      encode_dict[element1.GetString() + element2.GetString()] = element_encode;
      ++ct2;
    }
  }

  size_t num_triplets = type_size * type_size * type_size;
  size_t ct3 = 0;
  for (const auto &element1: type_set) {
    for (const auto &element2: type_set) {
      for (const auto &element3: type_set) {
        std::vector<double> element_encode(num_triplets, 0);
        element_encode[ct3] = 1.0;
        encode_dict[element1.GetString() + element2.GetString() + element3.GetString()] =
            element_encode;
        ++ct3;
      }
    }
  }
  return encode_dict;
}
inline std::vector<double> GetOneHotParametersFromMap(
    const std::vector<Element> & encode,
    const std::unordered_map<std::string, std::vector<double> > & one_hot_encode_hashmap,
    size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping){

  std::vector<double> res_encode;
  res_encode.reserve(354);

  for (const auto &cluster_vector: cluster_mapping) {
    std::vector<double> sum_of_list(
        static_cast<size_t>(std::pow(num_of_elements, cluster_vector[0].size())), 0);
    for (const auto &cluster: cluster_vector) {
      std::string cluster_type;
      for (auto index: cluster) {
        cluster_type += encode[index].GetString();
      }
      const auto &cluster_one_hot_encode = one_hot_encode_hashmap.at(cluster_type);
      std::transform(sum_of_list.begin(), sum_of_list.end(),
                     cluster_one_hot_encode.begin(),
                     sum_of_list.begin(),
                     std::plus<>());
    }
    auto cluster_vector_size = static_cast<double>( cluster_vector.size());
    std::for_each(sum_of_list.begin(),
                  sum_of_list.end(),
                  [cluster_vector_size](auto &n) { n /= cluster_vector_size; });

    std::move(sum_of_list.begin(), sum_of_list.end(), std::back_inserter(res_encode));
  }
  return res_encode;
}
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMMM(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMMM(
    const cfg::Config &config);

std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMM2(
    const cfg::Config &config);
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMM2(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);

std::vector<cfg::Lattice> GetSortedLatticeVectorPeriodic(
    const cfg::Config &config, size_t lattice_id);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingPeriodic(
    const cfg::Config &config);
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_LATTICECLUSTER_HPP_

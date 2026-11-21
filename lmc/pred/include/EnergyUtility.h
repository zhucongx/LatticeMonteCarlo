#ifndef LMC_LMC_PRED_INCLUDE_ENERGYUTILITY_H_
#define LMC_LMC_PRED_INCLUDE_ENERGYUTILITY_H_
#include "Config.h"
#include "ElementCluster.hpp"
#include "LatticeCluster.hpp"

#include <Eigen/Dense>
#include <array>
#include <boost/functional/hash.hpp>
#include <cstdint>
#include <limits>
#include <nlohmann/json.hpp>
#include <set>
#include <string>

namespace pred {
// clang-format off
inline Eigen::VectorXd JsonToEigenVector(const nlohmann::json &j) {
  const auto vec = j.get<std::vector<double>>();
  return Eigen::Map<const Eigen::VectorXd>(vec.data(), static_cast<Eigen::Index>(vec.size()));
}

inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
JsonToEigenMatrix(const nlohmann::json &j) {
  const auto vv = j.get<std::vector<std::vector<double>>>();
  if (vv.empty()) {
    return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>();
  }
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat(
      static_cast<Eigen::Index>(vv.size()), static_cast<Eigen::Index>(vv[0].size()));
  for (size_t row = 0; row < vv.size(); ++row) {
    mat.row(static_cast<Eigen::Index>(row)) = Eigen::Map<const Eigen::VectorXd>(
        vv[row].data(), static_cast<Eigen::Index>(vv[row].size()));
  }
  return mat;
}
// clang-format on
using Singlet_MMM_t = cfg::LatticeClusterMMM<1>;
using Pair_MMM_t = cfg::LatticeClusterMMM<2>;
using Singlet_MM2_t = cfg::LatticeClusterMM2<1>;
using Pair_MM2_t = cfg::LatticeClusterMM2<2>;

using Singlet_State_t = cfg::LatticeCluster<1>;
using Pair_State_t = cfg::LatticeCluster<2>;
using Triplet_State_t = cfg::LatticeCluster<3>;

std::unordered_map<std::string, std::vector<double>> GetOneHotEncodeHashmap(const std::set<Element> &element_set);
// Returns forward and backward sorted lattice lists
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMMM(const cfg::Config &config,
                                                                 const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMM2(const cfg::Config &config,
                                                                 const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t>>> GetAverageClusterParametersMappingMM2(const cfg::Config &config);
std::vector<std::vector<std::vector<size_t>>> GetAverageClusterParametersMappingMMM(const cfg::Config &config);

std::unordered_map<cfg::ElementCluster, size_t, boost::hash<cfg::ElementCluster>>
InitializeClusterHashMap(const std::set<Element> &element_set);
int GetLabel(const std::vector<size_t> &lattice_index_list, const cfg::Config &config);
// Returns forward and backward sorted lattice lists
std::vector<cfg::Lattice> GetSortedLatticeVectorStateOfPair(const cfg::Config &config,
                                                            const std::pair<size_t, size_t> &lattice_id_pair);
std::vector<cfg::Lattice> GetSortedLatticeVectorStateOfSite(const cfg::Config &config, size_t lattice_id);
std::vector<std::vector<std::vector<size_t>>> GetClusterParametersMappingStatePair(const cfg::Config &config);
std::vector<std::vector<std::vector<size_t>>> GetClusterParametersMappingStateSite(const cfg::Config &config);
std::vector<std::vector<std::vector<size_t>>>
GetClusterParametersMappingStatePairOf(const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_pair);
std::vector<std::vector<std::vector<size_t>>> GetClusterParametersMappingStateSiteOf(const cfg::Config &config,
                                                                                     size_t lattice_id);

void GetOneHotParametersFromMap(const std::vector<Element> &encode,
                                const std::unordered_map<std::string, std::vector<double>> &one_hot_encode_hashmap,
                                size_t num_of_elements,
                                const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping,
                                std::vector<double> &out_encode);
// Length of the vector GetOneHotParametersFromMap produces for this cluster mapping.
size_t GetOneHotEncodeSize(size_t num_of_elements,
                           const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping);
// Throws if an array read from the coefficient file does not have the size the model requires. Eigen checks no
// dimensions in release builds, so a mismatch would otherwise go on to produce meaningless energies.
void CheckModelSize(const std::string &predictor_filename,
                    const std::string &name,
                    size_t actual_size,
                    size_t expected_size);

struct ParametersQuartic {
  Eigen::VectorXd mu_x_mmm{};
  Eigen::VectorXd mu_x_mm2{};
  Eigen::VectorXd sigma_x_mmm{};
  Eigen::VectorXd sigma_x_mm2{};
  Eigen::MatrixXd U_mmm{};
  Eigen::MatrixXd U_mm2{};
  Eigen::VectorXd theta_D{};
  Eigen::VectorXd theta_Ks{};
  double mu_y_D{};
  double mu_y_Ks{};
  double sigma_y_D{};
  double sigma_y_Ks{};
};

struct ParametersE0 {
  Eigen::VectorXd mu_x_mmm{};
  Eigen::VectorXd sigma_x_mmm{};
  Eigen::MatrixXd U_mmm{};
  Eigen::VectorXd theta_e0{};
  double mu_y_e0{};
  double sigma_y_e0{};
};

struct ParametersDE {
  Eigen::VectorXd mu_x{};
  Eigen::VectorXd sigma_x{};
  Eigen::MatrixXd U{};

  Eigen::VectorXd theta{};
  double mu_y{};
  double sigma_y{};
};

// The clusters of a mapping as rows of one contiguous array, so that the loops which run on every evaluation read
// neighbouring memory instead of following a pointer per cluster.
struct FlatClusterMapping {
  struct Group {
    size_t label{};           // position of the group in the mapping, which is the label of its clusters
    size_t cluster_size{};    // number of sites of each of its clusters
    size_t num_clusters{};
    size_t offset{};          // of its first site index in site_indices
  };
  static constexpr size_t kMaxClusterSize = 3;
  std::vector<Group> groups{};
  std::vector<uint8_t> site_indices{};
};
FlatClusterMapping FlattenClusterMapping(const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping);

// Helper to map clusters to contiguous indices and store per-cluster total bond counts.
class ClusterIndexer {
 public:
  ClusterIndexer() = default;
  ClusterIndexer(const std::map<cfg::ElementCluster, int> &ordered_clusters,
                 std::vector<double> cluster_total_bonds);
  [[nodiscard]] size_t Size() const { return cluster_index_to_cluster_.size(); }
  [[nodiscard]] size_t GetIndex(const cfg::ElementCluster &cluster) const;
  // Small number that stands for the element in the GetIndex() below. Elements that are part of no cluster share one.
  [[nodiscard]] uint8_t GetElementCode(const Element element) const {
    return element_codes_[static_cast<uint8_t>(hash_value(element))];
  }
  // Index of the cluster with this label whose elements have the given codes, in any order. The result is that of
  // GetIndex(cfg::ElementCluster(label, elements)), but it costs a table lookup instead of allocating, sorting and
  // hashing a cluster. The energy models ask for several hundred clusters per evaluation.
  [[nodiscard]] size_t GetIndex(const size_t label, const uint8_t *codes, const size_t size) const {
    if (label < label_sizes_.size() && label_sizes_[label] == size) {
      size_t combined_code = 0;
      for (size_t i = 0; i < size; ++i) {
        combined_code = combined_code * num_codes_ + codes[i];
      }
      const auto index = lookup_[label * codes_per_label_ + combined_code];
      if (index != kNotFound) {
        return index;
      }
    }
    ThrowClusterNotFound();
  }
  [[nodiscard]] const std::vector<cfg::ElementCluster> &GetClusters() const { return cluster_index_to_cluster_; }
  [[nodiscard]] const std::vector<double> &GetTotalBonds() const { return cluster_total_bonds_; }
 private:
  static constexpr uint32_t kNotFound = std::numeric_limits<uint32_t>::max();
  [[noreturn]] static void ThrowClusterNotFound();
  void BuildLookup();

  std::vector<cfg::ElementCluster> cluster_index_to_cluster_{};
  std::vector<double> cluster_total_bonds_{};
  std::unordered_map<cfg::ElementCluster, size_t, boost::hash<cfg::ElementCluster>> cluster_to_index_{};
  // Lookup by codes: lookup_[label * codes_per_label_ + combined code] is the index of the cluster or kNotFound.
  size_t num_codes_{};
  size_t codes_per_label_{};
  std::vector<size_t> label_sizes_{};
  std::array<uint8_t, 256> element_codes_{};
  std::vector<uint32_t> lookup_{};
};

// For every cluster of the mapping, takes one from the count of the cluster its sites form with the start codes and
// adds one to the count of the cluster they form with the end codes. The codes are ClusterIndexer::GetElementCode() of
// the elements on the sites the mapping refers to, and `count_changes` has an entry per cluster index.
inline void AddClusterCountChanges(const FlatClusterMapping &mapping,
                                   const ClusterIndexer &indexer,
                                   const std::vector<uint8_t> &start_codes,
                                   const std::vector<uint8_t> &end_codes,
                                   std::vector<int> &count_changes) {
  std::array<uint8_t, FlatClusterMapping::kMaxClusterSize> cluster_start{};
  std::array<uint8_t, FlatClusterMapping::kMaxClusterSize> cluster_end{};
  for (const auto &group: mapping.groups) {
    const uint8_t *site_indices = mapping.site_indices.data() + group.offset;
    for (size_t n = 0; n < group.num_clusters; ++n, site_indices += group.cluster_size) {
      for (size_t i = 0; i < group.cluster_size; ++i) {
        cluster_start[i] = start_codes[site_indices[i]];
        cluster_end[i] = end_codes[site_indices[i]];
      }
      --count_changes[indexer.GetIndex(group.label, cluster_start.data(), group.cluster_size)];
      ++count_changes[indexer.GetIndex(group.label, cluster_end.data(), group.cluster_size)];
    }
  }
}

// Reusable thread-local buffers to avoid hot-path allocations.
// Each function provides a unique buffer for a specific purpose to prevent aliasing.

// Integer buffers for start/end state counting.
inline std::vector<int> &GetThreadLocalStartCountsBuffer() {
  thread_local std::vector<int> buffer;
  return buffer;
}
inline std::vector<int> &GetThreadLocalEndCountsBuffer() {
  thread_local std::vector<int> buffer;
  return buffer;
}

// Element codes of the sites of a neighbourhood before and after a change.
inline std::vector<uint8_t> &GetThreadLocalStartCodesBuffer() {
  thread_local std::vector<uint8_t> buffer;
  return buffer;
}
inline std::vector<uint8_t> &GetThreadLocalEndCodesBuffer() {
  thread_local std::vector<uint8_t> buffer;
  return buffer;
}

// Element buffers for start/end state cluster construction.
inline std::vector<Element> &GetThreadLocalElementStartBuffer() {
  thread_local std::vector<Element> buffer;
  return buffer;
}
inline std::vector<Element> &GetThreadLocalElementEndBuffer() {
  thread_local std::vector<Element> buffer;
  return buffer;
}

// Double-precision buffers for various encoding vectors.
inline std::vector<double> &GetThreadLocalDeEncodeBuffer() {
  thread_local std::vector<double> buffer;
  return buffer;
}
inline std::vector<double> &GetThreadLocalEncodeMMMBuffer() {
  thread_local std::vector<double> buffer;
  return buffer;
}
inline std::vector<double> &GetThreadLocalEncodeMM2ForwardBuffer() {
  thread_local std::vector<double> buffer;
  return buffer;
}
inline std::vector<double> &GetThreadLocalEncodeMM2BackwardBuffer() {
  thread_local std::vector<double> buffer;
  return buffer;
}
inline std::vector<double> &GetThreadLocalSumOfListBuffer() {
  thread_local std::vector<double> buffer;
  return buffer;
}

}    // namespace pred

#endif    //LMC_LMC_PRED_INCLUDE_ENERGYUTILITY_H_

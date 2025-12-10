
#ifndef LMC_LMC_PRED_INCLUDE_ENERGYUTILITY_H_
#define LMC_LMC_PRED_INCLUDE_ENERGYUTILITY_H_
#include "Config.h"
#include "ElementCluster.hpp"
#include "LatticeCluster.hpp"

#include <Eigen/Dense>
#include <boost/functional/hash.hpp>
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
std::vector<double>
GetOneHotParametersFromMap(const std::vector<Element> &encode,
                           const std::unordered_map<std::string, std::vector<double>> &one_hot_encode_hashmap,
                           size_t num_of_elements,
                           const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping);

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

// Helper to map clusters to contiguous indices and store per-cluster total bond counts.
class ClusterIndexer {
 public:
  ClusterIndexer() = default;
  ClusterIndexer(const std::map<cfg::ElementCluster, int> &ordered_clusters,
                 std::vector<double> cluster_total_bonds);
  [[nodiscard]] size_t Size() const { return cluster_index_to_cluster_.size(); }
  [[nodiscard]] size_t GetIndex(const cfg::ElementCluster &cluster) const;
  [[nodiscard]] const std::vector<cfg::ElementCluster> &GetClusters() const { return cluster_index_to_cluster_; }
  [[nodiscard]] const std::vector<double> &GetTotalBonds() const { return cluster_total_bonds_; }
 private:
  std::vector<cfg::ElementCluster> cluster_index_to_cluster_{};
  std::vector<double> cluster_total_bonds_{};
  std::unordered_map<cfg::ElementCluster, size_t, boost::hash<cfg::ElementCluster>> cluster_to_index_{};
};

// Reusable thread-local buffers to avoid hot-path allocations.
// Provides a thread-local primary buffer for integer vectors.
inline std::vector<int> &GetThreadLocalIntPrimaryBuffer() {
  thread_local std::vector<int> buffer;
  return buffer;
}
// Provides a thread-local secondary buffer for integer vectors (distinct from primary).
inline std::vector<int> &GetThreadLocalIntSecondaryBuffer() {
  thread_local std::vector<int> buffer;
  return buffer;
}
// Provides a thread-local buffer for double vectors.
inline std::vector<double> &GetThreadLocalDoubleBuffer() {
  thread_local std::vector<double> buffer;
  return buffer;
}
// Provides a thread-local primary buffer for Element vectors.
inline std::vector<Element> &GetThreadLocalElementPrimaryBuffer() {
  thread_local std::vector<Element> buffer;
  return buffer;
}
// Provides a thread-local secondary buffer for Element vectors (distinct from primary).
inline std::vector<Element> &GetThreadLocalElementSecondaryBuffer() {
  thread_local std::vector<Element> buffer;
  return buffer;
}
}    // namespace pred

#endif    //LMC_LMC_PRED_INCLUDE_ENERGYUTILITY_H_

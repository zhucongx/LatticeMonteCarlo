#include "VacancyMigrationPredictorE0.h"

#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <utility>

namespace pred {
VacancyMigrationPredictorE0::VacancyMigrationPredictorE0(const std::string &predictor_filename,
                                                               const cfg::Config &reference_config,
                                                               std::set<Element> element_set)
    : element_set_(std::move(element_set)),
      one_hot_encode_hash_map_(GetOneHotEncodeHashmap(element_set_)),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_state_(GetClusterParametersMappingStatePair(reference_config)) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(element_set_copy);
  ordered_map_ =
      std::map<cfg::ElementCluster, int>(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  if (!ifs) {
    throw std::runtime_error("Cannot open " + predictor_filename);
  }
  nlohmann::json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = JsonToEigenVector(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[Element(element)] = ParametersE0{
        JsonToEigenVector(parameters.at("mu_x_mmm")),
        JsonToEigenVector(parameters.at("sigma_x_mmm")),
        JsonToEigenMatrix(parameters.at("U_mmm")),
        JsonToEigenVector(parameters.at("theta_e0")),
        parameters.at("mu_e0"),
        parameters.at("sigma_e0")};
  }
#pragma omp parallel for default(none) shared(reference_config)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector = GetSortedLatticeVectorStateOfPair(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_state;
      std::transform(sorted_lattice_vector.begin(),
                     sorted_lattice_vector.end(),
                     std::back_inserter(lattice_id_vector_state),
                     [](const auto &lattice) {
                       return lattice.GetId();
                     });
#pragma omp critical
      { site_bond_cluster_state_hashmap_[{i, j}] = lattice_id_vector_state; }
      auto sorted_lattice_vector_mmm = GetSymmetricallySortedLatticeVectorMMM(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mmm;
      std::transform(sorted_lattice_vector_mmm.begin(),
                     sorted_lattice_vector_mmm.end(),
                     std::back_inserter(lattice_id_vector_mmm),
                     [](const auto &lattice) {
                       return lattice.GetId();
                     });
#pragma omp critical
      { site_bond_cluster_mmm_hashmap_[{i, j}] = lattice_id_vector_mmm; }
    }
  }
}

VacancyMigrationPredictorE0::~VacancyMigrationPredictorE0() = default;

std::pair<double, double> VacancyMigrationPredictorE0::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetBarrierAndDiffFromLatticeIdPair(config,
                                            {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
                                             config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}

double VacancyMigrationPredictorE0::GetDe(const cfg::Config &config,
                                             const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto &lattice_id_vector = site_bond_cluster_state_hashmap_.at(lattice_id_jump_pair);
  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);
  // #pragma omp parallel for default(none) shared(config, lattice_id_jump_pair, lattice_id_vector, migration_element, start_hashmap, end_hashmap)
  for (size_t label = 0; label < mapping_state_.size(); ++label) {
    const auto &cluster_vector = mapping_state_.at(label);
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto index: cluster) {
        size_t lattice_id = lattice_id_vector[index];
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id));
        if (lattice_id == lattice_id_jump_pair.first) {
          element_vector_end.push_back(migration_element);
          continue;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          element_vector_end.emplace_back(ElementName::X);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id));
      }
      auto cluster_start = cfg::ElementCluster(static_cast<int>(label), element_vector_start);
      auto cluster_end = cfg::ElementCluster(static_cast<int>(label), element_vector_end);
      // #pragma omp critical
      {
        start_hashmap[cluster_start]++;
        end_hashmap[cluster_end]++;
      }
    }
  }
  std::map<cfg::ElementCluster, int> ordered(ordered_map_);
  std::vector<double> de_encode;
  de_encode.reserve(ordered.size());
  static const std::vector<double> cluster_counter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
  // not necessary to parallelize this loop
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto start = static_cast<double>(start_hashmap.at(cluster));
    auto end = static_cast<double>(end_hashmap.at(cluster));
    auto total_bond = cluster_counter[static_cast<size_t>(cluster.GetLabel())];
    de_encode.push_back((end - start) / total_bond);
  }

  const Eigen::Map<const Eigen::VectorXd> encode_vec(de_encode.data(), static_cast<Eigen::Index>(de_encode.size()));
  return base_theta_.dot(encode_vec);
}

double VacancyMigrationPredictorE0::GetE0(const cfg::Config &config,
                                             const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto lattice_id_vector_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mmm.size());
  for (auto index: lattice_id_vector_mmm) {
    ele_vector.push_back(config.GetElementAtLatticeId(index));
  }
  auto &encode_mmm = GetThreadLocalEncodeMMMBuffer();
  GetOneHotParametersFromMap(ele_vector, one_hot_encode_hash_map_, element_set_.size(), mapping_mmm_, encode_mmm);
  const auto &[mu_x_mmm, sigma_x_mmm, U_mmm, theta_e0, mu_y_e0, sigma_y_e0] = element_parameters_hashmap_.at(migration_element);

  // Vectorized normalization and matrix operation using Eigen.
  Eigen::Map<Eigen::VectorXd> encode_vec(encode_mmm.data(), static_cast<Eigen::Index>(encode_mmm.size()));
  encode_vec -= mu_x_mmm;
  encode_vec = encode_vec.cwiseQuotient(sigma_x_mmm);

  const auto &U_mat = U_mmm;
  double logD = (U_mat * encode_vec).dot(theta_e0);
  logD *= sigma_y_e0;
  logD += mu_y_e0;
  return std::exp(logD);
}

std::pair<double, double> VacancyMigrationPredictorE0::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  double dE = GetDe(config, lattice_id_jump_pair);
  double e0 = GetE0(config, lattice_id_jump_pair);
  const auto Ea = std::max(0.0, e0 + dE / 2);
  return {Ea, dE};
}
}    // namespace pred

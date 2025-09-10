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
      auto base_theta_json = parameters.at("theta");
      base_theta_ = {};
      for (const auto &theta: base_theta_json) {
        base_theta_.emplace_back(theta.get<double>());
      }
      // base_theta_ = std::vector<double>(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[Element(element)] = ParametersE0{parameters.at("mu_x_mmm"),
                                                                 parameters.at("sigma_x_mmm"),
                                                                 parameters.at("U_mmm"),
                                                                 parameters.at("theta_e0"),
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

  double dE = 0;
  const size_t cluster_size = base_theta_.size();
  // not necessary to parallelize this loop
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += base_theta_[i] * de_encode[i];
  }
  return dE;
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
  auto encode_mmm = GetOneHotParametersFromMap(ele_vector, one_hot_encode_hash_map_, element_set_.size(), mapping_mmm_);
  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);
  const size_t old_size = element_parameters.mu_x_mmm.size();
  const size_t new_size = element_parameters.theta_e0.size();
  // not necessary to parallelize this loop
  for (size_t i = 0; i < old_size; ++i) {
    encode_mmm.at(i) -= element_parameters.mu_x_mmm.at(i);
    encode_mmm.at(i) /= element_parameters.sigma_x_mmm.at(i);
  }
  double logD = 0;
  // #pragma omp parallel for default(none) shared(encode_mmm, new_size, old_size, element_parameters) reduction(+:logD)
  // TODO(perf): Vectorize with Eigen or BLAS (gemv/noalias) to speed up PCA dot-products.
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mmm.at(i) * element_parameters.U_mmm.at(j).at(i);
    }
    auto logD_it = pca_dot * element_parameters.theta_e0.at(j);
    logD += logD_it;
  }
  logD *= element_parameters.sigma_y_e0;
  logD += element_parameters.mu_y_e0;
  return std::exp(logD);
}

std::pair<double, double> VacancyMigrationPredictorE0::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  double dE, e0;
  dE = GetDe(config, lattice_id_jump_pair);
  e0 = GetE0(config, lattice_id_jump_pair);
  const auto Ea = std::max(0.0, e0 + dE / 2);
  return {Ea, dE};
}
}    // namespace pred

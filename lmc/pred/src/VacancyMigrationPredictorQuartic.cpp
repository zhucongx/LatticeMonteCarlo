#include "VacancyMigrationPredictorQuartic.h"
#include <Eigen/Dense>
#include <utility>
#include <boost/range/combine.hpp>
#include <omp.h>
#include <nlohmann/json.hpp>

namespace pred {
VacancyMigrationPredictorQuartic::VacancyMigrationPredictorQuartic(const std::string &predictor_filename,
                                                                   const cfg::Config &reference_config,
                                                                   std::set<Element> element_set)
    : element_set_(std::move(element_set)),
      one_hot_encode_hash_map_(GetOneHotEncodeHashmap(element_set_)),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_mm2_(GetAverageClusterParametersMappingMM2(reference_config)),
      mapping_state_(GetClusterParametersMappingStatePair(reference_config)) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(element_set_copy);
  ordered_map_ = std::map<cfg::ElementCluster, int>(initialized_cluster_hashmap_.begin(),
                                                    initialized_cluster_hashmap_.end());

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  if (!ifs) {
    throw std::runtime_error("Cannot open " + predictor_filename);
  }
  nlohmann::json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters] : all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = JsonToEigenVector(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[Element(element)] = ParametersQuartic{
        JsonToEigenVector(parameters.at("mu_x_mmm")),
        JsonToEigenVector(parameters.at("mu_x_mm2")),
        JsonToEigenVector(parameters.at("sigma_x_mmm")),
        JsonToEigenVector(parameters.at("sigma_x_mm2")),
        JsonToEigenMatrix(parameters.at("U_mmm")),
        JsonToEigenMatrix(parameters.at("U_mm2")),
        JsonToEigenVector(parameters.at("theta_D")),
        JsonToEigenVector(parameters.at("theta_Ks")),
        parameters.at("mu_D"),
        parameters.at("mu_Ks"),
        parameters.at("sigma_D"),
        parameters.at("sigma_Ks")
    };
  }
#pragma omp parallel for default(none) shared(reference_config)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j : reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector =
          GetSortedLatticeVectorStateOfPair(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_state;
      std::transform(sorted_lattice_vector.begin(), sorted_lattice_vector.end(),
                     std::back_inserter(lattice_id_vector_state),
                     [](const auto &lattice) { return lattice.GetId(); });
#pragma omp critical
      {
        site_bond_cluster_state_hashmap_[{i, j}] = lattice_id_vector_state;
      }
      auto sorted_lattice_vector_mmm =
          GetSymmetricallySortedLatticeVectorMMM(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mmm;
      std::transform(sorted_lattice_vector_mmm.begin(), sorted_lattice_vector_mmm.end(),
                     std::back_inserter(lattice_id_vector_mmm),
                     [](const auto &lattice) { return lattice.GetId(); });
#pragma omp critical
      {
        site_bond_cluster_mmm_hashmap_[{i, j}] = lattice_id_vector_mmm;
      }
      auto sorted_lattice_vector_mm2 =
          GetSymmetricallySortedLatticeVectorMM2(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mm2;
      std::transform(sorted_lattice_vector_mm2.begin(), sorted_lattice_vector_mm2.end(),
                     std::back_inserter(lattice_id_vector_mm2),
                     [](const auto &lattice) { return lattice.GetId(); });
#pragma omp critical
      {
        site_bond_cluster_mm2_hashmap_[{i, j}] = lattice_id_vector_mm2;
      }
      // TODO(arch): Consider storing these neighbor lists in contiguous arrays with
      // deterministic indices (site*12 + nn_idx). This avoids unordered_map<pair>
      // lookups on hot paths and reduces bounds-checked .at() overhead.
    }
  }
}
VacancyMigrationPredictorQuartic::~VacancyMigrationPredictorQuartic() = default;

std::pair<double, double> VacancyMigrationPredictorQuartic::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetBarrierAndDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
double VacancyMigrationPredictorQuartic::GetDe(const cfg::Config &config,
                                               const std::pair<size_t,
                                                               size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto &lattice_id_vector = site_bond_cluster_state_hashmap_.at(lattice_id_jump_pair);
  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);
  // TODO(perf): These per-call maps/vectors allocate each time. Replace with
  // preallocated thread_local fixed arrays to eliminate allocation/tree overhead.
// #pragma omp parallel for default(none) shared(config, lattice_id_jump_pair, lattice_id_vector, migration_element, start_hashmap, end_hashmap)
  for (size_t label = 0; label < mapping_state_.size(); ++label) {
    const auto &cluster_vector = mapping_state_.at(label);
    for (const auto &cluster : cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto index : cluster) {
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
  static const std::vector<double>
      cluster_counter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
  // not necessary to parallelize this loop
  for (const auto &cluster_count : ordered) {
    const auto &cluster = cluster_count.first;
    auto start = static_cast<double>(start_hashmap.at(cluster));
    auto end = static_cast<double>(end_hashmap.at(cluster));
    auto total_bond = cluster_counter[static_cast<size_t>(cluster.GetLabel())];
    de_encode.push_back((end - start) / total_bond);
  }

  const Eigen::Map<const Eigen::VectorXd> encode_vec(de_encode.data(), static_cast<Eigen::Index>(de_encode.size()));
  return base_theta_.dot(encode_vec);
}
double VacancyMigrationPredictorQuartic::GetKs(const cfg::Config &config,
                                               const std::pair<size_t,
                                                               size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);

  auto lattice_id_vector_mm2_forward =
      site_bond_cluster_mm2_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector_forward{};
  ele_vector_forward.reserve(lattice_id_vector_mm2_forward.size());
  // TODO(perf): Convert these per-call vectors to preallocated thread_local arrays
  // to avoid allocations inside the predictor hot path.
  for (auto index : lattice_id_vector_mm2_forward) {
    ele_vector_forward.push_back(config.GetElementAtLatticeId(index));
  }
  auto encode_mm2_forward = GetOneHotParametersFromMap(ele_vector_forward,
                                                       one_hot_encode_hash_map_,
                                                       element_set_.size(),
                                                       mapping_mm2_);

  auto lattice_id_vector_mm2_backward =
      site_bond_cluster_mm2_hashmap_.at({lattice_id_jump_pair.second, lattice_id_jump_pair.first});
  std::vector<Element> ele_vector_backward{};
  ele_vector_backward.reserve(lattice_id_vector_mm2_backward.size());
  for (auto index : lattice_id_vector_mm2_backward) {
    ele_vector_backward.push_back(config.GetElementAtLatticeId(index));
  }
  const auto encode_mm2_backward = GetOneHotParametersFromMap(ele_vector_backward,
                                                        one_hot_encode_hash_map_,
                                                        element_set_.size(),
                                                        mapping_mm2_);

  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);
  // Vectorized normalization using Eigen.
  Eigen::Map<Eigen::VectorXd> encode_vec(encode_mm2_forward.data(), static_cast<Eigen::Index>(encode_mm2_forward.size()));
  const Eigen::Map<const Eigen::VectorXd> backward_vec(encode_mm2_backward.data(),  static_cast<Eigen::Index>(encode_mm2_backward.size()));
  encode_vec += backward_vec;
  encode_vec -= element_parameters.mu_x_mm2;
  encode_vec = encode_vec.cwiseQuotient(element_parameters.sigma_x_mm2);

  const auto &U_mat = element_parameters.U_mm2;
  double logKs = (U_mat * encode_vec).dot(element_parameters.theta_Ks);
  logKs *= element_parameters.sigma_y_Ks;
  logKs += element_parameters.mu_y_Ks;
  return std::exp(logKs);
}
double VacancyMigrationPredictorQuartic::GetD(const cfg::Config &config,
                                              const std::pair<size_t,
                                                              size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto lattice_id_vector_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mmm.size());
  // TODO(perf): Replace per-call vector alloc with thread_local fixed storage
  // to eliminate allocator overhead.
  for (auto index : lattice_id_vector_mmm) {
    ele_vector.push_back(config.GetElementAtLatticeId(index));
  }
  auto encode_mmm = GetOneHotParametersFromMap(ele_vector,
                                               one_hot_encode_hash_map_,
                                               element_set_.size(),
                                               mapping_mmm_);

  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);

  // Vectorized normalization using Eigen.
  Eigen::Map<Eigen::VectorXd> encode_vec(encode_mmm.data(), static_cast<Eigen::Index>(encode_mmm.size()));
  encode_vec -= element_parameters.mu_x_mmm;
  encode_vec = encode_vec.cwiseQuotient(element_parameters.sigma_x_mmm);

  const auto &U_mat = element_parameters.U_mmm;
  double logD = (U_mat * encode_vec).dot(element_parameters.theta_D);
  logD *= element_parameters.sigma_y_D;
  logD += element_parameters.mu_y_D;
  return std::exp(logD);
}
std::pair<double, double> VacancyMigrationPredictorQuartic::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  double dE, D, Ks;
#pragma omp parallel sections default(none) shared(config, lattice_id_jump_pair, dE, D, Ks)
  {
#pragma omp section
    {
      dE = GetDe(config, lattice_id_jump_pair);
    }
#pragma omp section
    {
      D = GetD(config, lattice_id_jump_pair);
    }
#pragma omp section
    {
      Ks = GetKs(config, lattice_id_jump_pair);
    }
  }
  // TODO(perf): When already inside an outer OMP region (e.g., FirstOmp), gate this
  // with omp_in_parallel() and run sequential to avoid nested-team overhead. Conversely,
  // if outer loop is sequential, keep these sections to parallelize dE/D/Ks.
  // const auto dE = GetDe(config, lattice_id_jump_pair);
  // const auto D = GetD(config, lattice_id_jump_pair);
  // const auto Ks = GetKs(config, lattice_id_jump_pair);
  const auto b = 4 * dE / (D * D * D);
  const auto a = Ks / (4 * D * D);
  const auto c = (9 * b * b - 16 * a * a * D * D) / (32 * a);
  const auto delta = std::sqrt(std::abs(9 * b * b - 32 * a * c));
  const auto Ea = (3 * b + delta) * (3 * b + delta) * (3 * b * b - 16 * a * c + b * delta)
      / std::pow(a, 3) / 2048;
  return {Ea, dE};
}
} // pred

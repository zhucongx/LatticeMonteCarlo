#include "VacancyMigrationPredictorQuartic.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
VacancyMigrationPredictorQuartic::VacancyMigrationPredictorQuartic(const std::string &predictor_filename,
                                                                   const cfg::Config &reference_config,
                                                                   std::set<Element> type_set)
    : type_set_(std::move(type_set)),
      one_hot_encode_hash_map_(GetOneHotEncodeHashmap(type_set_)),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_mm2_(GetAverageClusterParametersMappingMM2(reference_config)),
      mapping_state_(GetClusterParametersMappingState(reference_config)) {
  auto type_set_copy(type_set_);
  type_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(type_set_copy);

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = std::vector<double>(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[Element(element)] = ParametersQuartic{
        parameters.at("mu_x_mmm"),
        parameters.at("mu_x_mm2"),
        parameters.at("sigma_x_mmm"),
        parameters.at("sigma_x_mm2"),
        parameters.at("U_mmm"),
        parameters.at("U_mm2"),
        parameters.at("theta_D"),
        parameters.at("theta_Ks"),
        parameters.at("mu_D"),
        parameters.at("mu_Ks"),
        parameters.at("sigma_D"),
        parameters.at("sigma_Ks")
    };
  }

  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector =
          GetSortedLatticeVectorState(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_state;
      std::transform(sorted_lattice_vector.begin(), sorted_lattice_vector.end(),
                     std::back_inserter(lattice_id_vector_state),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_state_hashmap_[{i, j}] = lattice_id_vector_state;
      auto sorted_lattice_vector_mmm =
          GetSymmetricallySortedLatticeVectorMMM(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mmm;
      std::transform(sorted_lattice_vector_mmm.begin(), sorted_lattice_vector_mmm.end(),
                     std::back_inserter(lattice_id_vector_mmm),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_mmm_hashmap_[{i, j}] = lattice_id_vector_mmm;
      auto sorted_lattice_vector_mm2 =
          GetSymmetricallySortedLatticeVectorMM2(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mm2;
      std::transform(sorted_lattice_vector_mm2.begin(), sorted_lattice_vector_mm2.end(),
                     std::back_inserter(lattice_id_vector_mm2),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_mm2_hashmap_[{i, j}] = lattice_id_vector_mm2;
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

  int label = 0;
  for (const auto &cluster_vector: mapping_state_) {
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
      start_hashmap[cfg::ElementCluster(label, element_vector_start)]++;
      end_hashmap[cfg::ElementCluster(label, element_vector_end)]++;
    }
    label++;
  }

  std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  std::vector<double> de_encode;
  de_encode.reserve(ordered.size());
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto start = static_cast<double>(start_hashmap.at(cluster));
    auto end = static_cast<double>(end_hashmap.at(cluster));
    double total_bond{};
    switch (cluster.GetLabel()) {
      case 0:total_bond = 256;
        break;
      case 1: total_bond = 1536;
        break;
      case 2: total_bond = 768;
        break;
      case 3: total_bond = 3072;
        break;
      case 4: total_bond = 2048;
        break;
      case 5: total_bond = 3072;
        break;
      case 6: total_bond = 6144;
        break;
      case 7: total_bond = 6144;
        break;
      case 8: total_bond = 6144;
        break;
      case 9: total_bond = 6144;
        break;
      case 10: total_bond = 2048;
        break;
    }
    de_encode.push_back((end - start) / total_bond);
  }
  double dE = 0;
  const size_t cluster_size = base_theta_.size();
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += base_theta_[i] * de_encode[i];
  }
  return dE;
}
double VacancyMigrationPredictorQuartic::GetKs(const cfg::Config &config,
                                               const std::pair<size_t,
                                                               size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);

  auto lattice_id_vector_mm2_forward =
      site_bond_cluster_mm2_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector_forward{};
  ele_vector_forward.reserve(lattice_id_vector_mm2_forward.size());
  for (auto index: lattice_id_vector_mm2_forward) {
    ele_vector_forward.push_back(config.GetElementAtLatticeId(index));
  }
  auto encode_mm2_forward = GetOneHotParametersFromMap(ele_vector_forward,
                                                       one_hot_encode_hash_map_,
                                                       type_set_.size(),
                                                       mapping_mm2_);

  auto lattice_id_vector_mm2_backward =
      site_bond_cluster_mm2_hashmap_.at({lattice_id_jump_pair.second, lattice_id_jump_pair.first});
  std::vector<Element> ele_vector_backward{};
  ele_vector_backward.reserve(lattice_id_vector_mm2_backward.size());
  for (auto index: lattice_id_vector_mm2_backward) {
    ele_vector_backward.push_back(config.GetElementAtLatticeId(index));
  }
  auto encode_mm2_backward = GetOneHotParametersFromMap(ele_vector_backward,
                                                        one_hot_encode_hash_map_,
                                                        type_set_.size(),
                                                        mapping_mm2_);

  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);

  const auto &mu_x = element_parameters.mu_x_mm2;
  const auto &sigma_x = element_parameters.sigma_x_mm2;
  const auto &U = element_parameters.U_mm2;
  const auto &theta = element_parameters.theta_Ks;
  const auto mu_y = element_parameters.mu_y_Ks;
  const auto sigma_y = element_parameters.sigma_y_Ks;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode_mm2_forward.at(i) += encode_mm2_backward.at(i);
    encode_mm2_forward.at(i) -= mu_x.at(i);
    encode_mm2_forward.at(i) /= sigma_x.at(i);
  }
  double logKs = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mm2_forward.at(i) * U.at(j).at(i);
    }
    logKs += pca_dot * theta.at(j);
  }
  logKs *= sigma_y;
  logKs += mu_y;
  return std::exp(logKs);
}
double VacancyMigrationPredictorQuartic::GetD(const cfg::Config &config,
                                              const std::pair<size_t,
                                                              size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto lattice_id_vector_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mmm.size());
  for (auto index: lattice_id_vector_mmm) {
    ele_vector.push_back(config.GetElementAtLatticeId(index));
  }
  auto encode_mmm = GetOneHotParametersFromMap(ele_vector,
                                               one_hot_encode_hash_map_,
                                               type_set_.size(),
                                               mapping_mmm_);
  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);

  const auto &mu_x = element_parameters.mu_x_mmm;
  const auto &sigma_x = element_parameters.sigma_x_mmm;
  const auto &U = element_parameters.U_mmm;
  const auto &theta = element_parameters.theta_D;
  const auto mu_y = element_parameters.mu_y_D;
  const auto sigma_y = element_parameters.sigma_y_D;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode_mmm.at(i) -= mu_x.at(i);
    encode_mmm.at(i) /= sigma_x.at(i);
  }
  double logD = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mmm.at(i) * U.at(j).at(i);
    }
    logD += pca_dot * theta.at(j);
  }
  logD *= sigma_y;
  logD += mu_y;
  return std::exp(logD);
}
std::pair<double, double> VacancyMigrationPredictorQuartic::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto dE = GetDe(config, lattice_id_jump_pair);
  const auto D = GetD(config, lattice_id_jump_pair);
  const auto Ks = GetKs(config, lattice_id_jump_pair);
  const auto b = 4 * dE / (D * D * D);
  const auto a = Ks / (4 * D * D);
  const auto c = (9 * b * b - 16 * a * a * D * D) / (32 * a);
  const auto delta = std::sqrt(std::abs(9 * b * b - 32 * a * c));
  const auto Ea = (3 * b + delta) * (3 * b + delta) * (3 * b * b - 16 * a * c + b * delta)
      / std::pow(a, 3) / 2048;
  return {Ea, dE};
}
} // namespace pred

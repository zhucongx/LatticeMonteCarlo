#include "EnergyPredictorQuartic.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace preds {

EnergyPredictorQuartic::EnergyPredictorQuartic(const std::string &predictor_filename,
                                               const cfg::Config &reference_config,
                                               const std::set<Element> &type_set)
    : EnergyPredictor(type_set),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_mm2_(GetAverageClusterParametersMappingMM2(reference_config)) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    element_parameters_hashmap_[Element(element)] = ParametersQuartic{
        parameters.at("mu_x_mmm"),
        parameters.at("mu_x_mm2"),
        parameters.at("sigma_x_mmm"),
        parameters.at("sigma_x_mm2"),
        parameters.at("U_mmm"),
        parameters.at("U_mm2"),
        parameters.at("theta_a"),
        parameters.at("theta_b"),
        parameters.at("theta_c"),
        parameters.at("mu_a"),
        parameters.at("mu_b"),
        parameters.at("mu_c"),
        parameters.at("sigma_a"),
        parameters.at("sigma_b"),
        parameters.at("sigma_c")
    };
  }
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
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

EnergyPredictorQuartic::~EnergyPredictorQuartic() = default;

std::pair<double, double> EnergyPredictorQuartic::GetAAndC(const cfg::Config &config,
                                                           const std::pair<size_t,
                                                                           size_t> &lattice_id_jump_pair,
                                                           Element migration_element) const {
  auto lattice_id_vector_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mmm.size());
  for (auto index: lattice_id_vector_mmm) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementType::X) {
      ele_vector.push_back(migration_element);
      continue;
    }
    ele_vector.push_back(this_element);
  }
  auto encode_mmm = pred::GetOneHotParametersFromMap(ele_vector,
                                                     one_hot_encode_hash_map_,
                                                     type_set_.size(),
                                                     mapping_mmm_);
  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);

  const auto &mu_x_mmm = element_parameters.mu_x_mmm;
  const auto &sigma_x_mmm = element_parameters.sigma_x_mmm;
  const auto &U_mmm = element_parameters.U_mmm;
  const auto &theta_a = element_parameters.theta_a;
  const auto mu_a = element_parameters.mu_a;
  const auto sigma_a = element_parameters.sigma_a;
  const auto &theta_c = element_parameters.theta_c;
  const auto mu_c = element_parameters.mu_c;
  const auto sigma_c = element_parameters.sigma_c;

  const size_t old_size = mu_x_mmm.size();
  const size_t new_size = theta_a.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode_mmm[i] -= mu_x_mmm[i];
    encode_mmm[i] /= sigma_x_mmm[i];
  }

  double a = 0, c = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mmm[i] * U_mmm[j][i];
    }
    a += pca_dot * theta_a[j];
    c += pca_dot * theta_c[j];
  }
  a *= sigma_a;
  a += mu_a;

  c *= sigma_c;
  c += mu_c;
  return {a, c};
}
double EnergyPredictorQuartic::GetB(const cfg::Config &config,
                                    const std::pair<size_t,
                                                    size_t> &lattice_id_jump_pair,
                                    Element migration_element) const {
  auto lattice_id_vector_mm2 = site_bond_cluster_mm2_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mm2.size());
  for (auto index: lattice_id_vector_mm2) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementType::X) {
      ele_vector.push_back(migration_element);
      continue;
    }
    ele_vector.push_back(this_element);
  }
  auto encode_mm2 = pred::GetOneHotParametersFromMap(ele_vector,
                                                     one_hot_encode_hash_map_,
                                                     type_set_.size(),
                                                     mapping_mm2_);
  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);

  const auto &mu_x_mm2 = element_parameters.mu_x_mm2;
  const auto &sigma_x_mm2 = element_parameters.sigma_x_mm2;
  const auto &U_mm2 = element_parameters.U_mm2;
  const auto &theta_b = element_parameters.theta_b;
  const auto mu_b = element_parameters.mu_b;
  const auto sigma_b = element_parameters.sigma_b;

  const size_t old_size = mu_x_mm2.size();
  const size_t new_size = theta_b.size();
  for (size_t i = 0; i < old_size; ++i) {
    encode_mm2[i] -= mu_x_mm2[i];
    encode_mm2[i] /= sigma_x_mm2[i];
  }

  double b = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mm2[i] * U_mm2[j][i];
    }
    b += pca_dot * theta_b[j];
  }
  b *= sigma_b;
  b += mu_b;
  return b;
}
std::pair<double, double> EnergyPredictorQuartic::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto[a1, c1] = GetAAndC(config,
                          lattice_id_jump_pair,
                          migration_element);
  auto[a2, c2] = GetAAndC(config,
                          {lattice_id_jump_pair.second, lattice_id_jump_pair.first},
                          migration_element);
  auto b1 = GetB(config,
                 lattice_id_jump_pair,
                 migration_element);
  auto b2 = GetB(config,
                 {lattice_id_jump_pair.second, lattice_id_jump_pair.first},
                 migration_element);
  auto a = (a1 + a2) / 2;
  auto b = (b1 - b2) / 2;
  auto c = (c1 + c2) / 2;

  // auto delta = std::sqrt(9 * b * b - 32 * a * c);
  auto delta = std::sqrt(std::abs(9 * b * b - 32 * a * c));
  auto dE = b * std::pow((delta / a), 3) / 256;
  auto Ea = (3 * b + delta) * (3 * b + delta) * (3 * b * b - 16 * a * c + b * delta)
      / std::pow(a, 3) / 2048;
  return {Ea, dE};
}

} // namespace pred

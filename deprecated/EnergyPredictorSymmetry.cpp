#include "EnergyPredictorSymmetry.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace pred {

EnergyPredictorSymmetry::EnergyPredictorSymmetry(const std::string &predictor_filename,
                                                 const cfg::Config &reference_config,
                                                 const std::set<Element> &type_set)
    : one_hot_encode_hash_map_(GetOneHotEncodeHashmap(type_set)),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_mm2_(GetAverageClusterParametersMappingMM2(reference_config)),
      type_set_(type_set) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    e0_element_parameters_hashmap_[Element(element)] = ParametersE0{
        parameters.at("e0_mu_x"),
        parameters.at("e0_sigma_x"),
        parameters.at("e0_mu_y"),
        parameters.at("e0_sigma_y"),
        parameters.at("e0_U"),
        parameters.at("e0_theta"),
    };
    dE_element_parameters_hashmap_[Element(element)] = ParametersDE{
        parameters.at("dE_mu_x"),
        parameters.at("dE_sigma_x"),
        parameters.at("dE_mu_y"),
        parameters.at("dE_sigma_y"),
        parameters.at("dE_U"),
        parameters.at("dE_theta"),
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
EnergyPredictorSymmetry::~EnergyPredictorSymmetry() = default;
std::pair<double, double> EnergyPredictorSymmetry::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {

  return GetBarrierAndDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
std::pair<double, double> EnergyPredictorSymmetry::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto e0 = GetE0(config,
                  lattice_id_jump_pair,
                  migration_element);

  auto ef = GetE(config,
                 lattice_id_jump_pair,
                 migration_element);
  auto eb = GetE(config,
                 {lattice_id_jump_pair.second, lattice_id_jump_pair.first},
                 migration_element);
  auto dE = eb - ef;
  return {e0 + dE / 2, dE};
}
double EnergyPredictorSymmetry::GetE0(const cfg::Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                      Element migration_element) const {
  auto lattice_id_vector_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mmm.size());
  for (auto index: lattice_id_vector_mmm) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementName::X) {
      ele_vector.push_back(migration_element);
      continue;
    }
    ele_vector.push_back(this_element);
  }
  auto encode_mmm = GetOneHotParametersFromMap(ele_vector,
                                               one_hot_encode_hash_map_,
                                               type_set_.size(),
                                               mapping_mmm_);
  const auto &element_parameters = e0_element_parameters_hashmap_.at(migration_element);

  const auto &mu_x = element_parameters.mu_x;
  const auto &sigma_x = element_parameters.sigma_x;
  const auto mu_y = element_parameters.mu_y;
  const auto sigma_y = element_parameters.sigma_y;

  const auto &U = element_parameters.U;
  const auto &theta = element_parameters.theta;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode_mmm.at(i) -= mu_x.at(i);
    encode_mmm.at(i) /= sigma_x.at(i);
  }
  double e0 = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mmm.at(i) * U.at(j).at(i);
    }
    e0 += pca_dot * theta.at(j);
  }
  e0 *= sigma_y;
  e0 += mu_y;
  return std::exp(e0);
}
double EnergyPredictorSymmetry::GetE(const cfg::Config &config,
                                     const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                     Element migration_element) const {
  auto lattice_id_vector_mm2 = site_bond_cluster_mm2_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mm2.size());
  for (auto index: lattice_id_vector_mm2) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementName::X) {
      ele_vector.push_back(migration_element);
      continue;
    }
    ele_vector.push_back(this_element);
  }
  auto encode_mm2 = GetOneHotParametersFromMap(ele_vector,
                                               one_hot_encode_hash_map_,
                                               type_set_.size(),
                                               mapping_mm2_);
  const auto &element_parameters = dE_element_parameters_hashmap_.at(migration_element);

  const auto &mu_x = element_parameters.mu_x;
  const auto &sigma_x = element_parameters.sigma_x;
  const auto mu_y = element_parameters.mu_y;
  const auto sigma_y = element_parameters.sigma_y;

  const auto &U = element_parameters.U;
  const auto &theta = element_parameters.theta;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode_mm2.at(i) -= mu_x.at(i);
    encode_mm2.at(i) /= sigma_x.at(i);
  }
  double e = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mm2.at(i) * U.at(j).at(i);
    }
    e += pca_dot * theta.at(j);
  }
  e *= sigma_y;
  e += mu_y;
  return e;
}
} // namespace pred

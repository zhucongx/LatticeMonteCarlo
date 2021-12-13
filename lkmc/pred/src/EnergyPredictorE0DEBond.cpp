#include "EnergyPredictorE0DEBond.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
EnergyPredictorE0DEBond::EnergyPredictorE0DEBond(const std::string &predictor_filename,
                                                 const cfg::Config &reference_config,
                                                 const std::set<Element> &type_set)
    : EnergyPredictor(type_set),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      initialized_hashmap_(InitializeHashMap(type_set)) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Bond") {
      theta_bonds = std::vector<double>(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[Element(element)] = ParametersE0{
        parameters.at("mu_x"),
        parameters.at("sigma_x"),
        parameters.at("mu_y"),
        parameters.at("sigma_y"),
        parameters.at("U"),
        parameters.at("theta"),
    };

  }
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector_mmm =
          GetSymmetricallySortedLatticeVectorMMM(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector;
      std::transform(sorted_lattice_vector_mmm.begin(), sorted_lattice_vector_mmm.end(),
                     std::back_inserter(lattice_id_vector),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_mmm_hashmap_[{i, j}] = lattice_id_vector;
    }
  }
}

EnergyPredictorE0DEBond::~EnergyPredictorE0DEBond() = default;

double EnergyPredictorE0DEBond::GetE0(const cfg::Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair,
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
  auto encode = pred::GetOneHotParametersFromMap(ele_vector,
                                                 one_hot_encode_hash_map_,
                                                 type_set_.size(),
                                                 mapping_mmm_);
  const auto &element_parameters = element_parameters_hashmap_.at(migration_element);

  const auto &mu_x = element_parameters.mu_x;
  const auto &sigma_x = element_parameters.sigma_x;
  const auto mu_y = element_parameters.mu_y;
  const auto sigma_y = element_parameters.sigma_y;

  const auto &U = element_parameters.U;
  const auto &theta = element_parameters.theta;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode[i] -= mu_x[i];
    encode[i] /= sigma_x[i];
  }
  double e0 = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode[i] * U[j][i];
    }
    e0 += pca_dot * theta[j];
  }
  e0 *= sigma_y;
  e0 += mu_y;
  return e0;
}
double EnergyPredictorE0DEBond::GetDE(const cfg::Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto bond_change_vector = GetBondChange(config, lattice_id_jump_pair, initialized_hashmap_);
  const size_t bond_size = theta_bonds.size();
  double dE = 0;
  for (size_t i = 0; i < bond_size; ++i) {
    dE += theta_bonds[i] * bond_change_vector[i];
  }
  return dE;
}
std::pair<double, double> EnergyPredictorE0DEBond::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto e01 = GetE0(config,
                   lattice_id_jump_pair,
                   migration_element);
  auto e02 = GetE0(config,
                   {lattice_id_jump_pair.second, lattice_id_jump_pair.first},
                   migration_element);
  auto e0 = (e01 + e02) / 2;

  auto dE = GetDE(config, lattice_id_jump_pair);

  return {e0 + dE / 2, dE};
}

} // namespace pred
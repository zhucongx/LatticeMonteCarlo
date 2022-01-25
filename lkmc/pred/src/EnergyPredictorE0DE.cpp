#include "EnergyPredictorE0DE.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {



EnergyPredictorE0DE::EnergyPredictorE0DE(const std::string &predictor_filename,
                                         const cfg::Config &reference_config,
                                         const std::set<Element> &type_set)
    : EnergyPredictor(type_set),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_periodic_(GetAverageClusterParametersMappingPeriodic(reference_config)) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "dE") {
      parameters_dE_ = ParametersDE{
          parameters.at("mu_x"),
          parameters.at("sigma_x"),
          parameters.at("mu_y"),
          parameters.at("sigma_y"),
          parameters.at("theta"),
      };
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
    auto sorted_lattice_vector_periodic = GetSortedLatticeVectorPeriodic(reference_config, i);
    std::vector<size_t> lattice_id_vector;
    std::transform(sorted_lattice_vector_periodic.begin(), sorted_lattice_vector_periodic.end(),
                   std::back_inserter(lattice_id_vector),
                   [](const auto &lattice) { return lattice.GetId(); });
    site_cluster_periodic_hashmap_[i] = lattice_id_vector;
  }
}

EnergyPredictorE0DE::~EnergyPredictorE0DE() = default;

double EnergyPredictorE0DE::GetE0(const cfg::Config &config,
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
double EnergyPredictorE0DE::GetE(const cfg::Config &config,
                                 size_t lattice_id,
                                 Element migration_element) const {
  auto lattice_id_vector_periodic = site_cluster_periodic_hashmap_.at(lattice_id);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_periodic.size());
  for (auto index: lattice_id_vector_periodic) {
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
                                                 mapping_periodic_);

  const auto &mu_x = parameters_dE_.mu_x;
  const auto &sigma_x = parameters_dE_.sigma_x;
  const auto mu_y = parameters_dE_.mu_y;
  const auto sigma_y = parameters_dE_.sigma_y;

  const auto &theta = parameters_dE_.theta;

  const size_t size = theta.size();
  for (size_t i = 0; i < size; ++i) {
    encode[i] -= mu_x[i];
    encode[i] /= sigma_x[i];
  }

  double e = 0;
  for (size_t i = 0; i < size; ++i) {
    e += theta[i] * encode[i];
  }
  e *= sigma_y;
  e += mu_y;
  return e;
}
std::pair<double, double> EnergyPredictorE0DE::GetBarrierAndDiffFromLatticeIdPair(
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

  auto dE = GetE(config, lattice_id_jump_pair.second, migration_element)
      - GetE(config, lattice_id_jump_pair.first, migration_element);

  return {e0 + dE / 2, dE};
}

} // namespace pred
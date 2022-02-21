#include "EnergyPredictorE0DEState.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {

EnergyPredictorE0DEState::EnergyPredictorE0DEState(const std::string &predictor_filename,
                                                   const cfg::Config &reference_config,
                                                   const std::set<Element> &type_set) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters]: all_parameters.items()) {
    element_theta_[Element(element)] = std::vector<double>(parameters.at("theta"));
  }
  for (const auto &element: type_set) {
    auto type_set_copy(type_set);
    type_set_copy.emplace(ElementType::X);
    type_set_copy.insert(element.GetPseudo());
    element_initialized_cluster_hashmap_[element] = InitializeClusterHashMap(type_set_copy);
  }

  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector =
          GetSortedLatticeVectorState(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector;
      std::transform(sorted_lattice_vector.begin(), sorted_lattice_vector.end(),
                     std::back_inserter(lattice_id_vector),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_mapping_hashmap_[{i, j}] =
          GetClusterParametersMappingState(reference_config, {i, j});
    }
  }
}
EnergyPredictorE0DEState::~EnergyPredictorE0DEState() = default;
std::pair<double, double> EnergyPredictorE0DEState::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto &cluster_mapping = site_bond_mapping_hashmap_.at(lattice_id_jump_pair);
  const auto &initialized_cluster_hashmap
      = element_initialized_cluster_hashmap_.at(migration_element);

  auto[de_encode, e0_encode] = GetEncodesFromMapState(config,
                                                      lattice_id_jump_pair,
                                                      initialized_cluster_hashmap,
                                                      cluster_mapping);
  const auto &theta = element_theta_.at(migration_element);
  double dE = 0, e0 = 0;

  const size_t cluster_size = theta.size();
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += theta[i] * de_encode[i];
    e0 += theta[i] * e0_encode[i];
  }
  return {e0 + dE / 2, dE};
}
} // namespace pred
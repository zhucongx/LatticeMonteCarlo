#include "EnergyPredictor.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <boost/functional/hash.hpp>

#include <omp.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace pred {

EnergyPredictor::EnergyPredictor(const std::string &predictor_filename,
                                 const cfg::Config &reference_config,
                                 const std::set<Element> &type_set) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = std::vector<double>(parameters.at("theta"));
      continue;
    }
    element_theta_[Element(element)] = std::vector<double>(parameters.at("theta"));
  }
  for (const auto &element: type_set) {
    auto type_set_copy(type_set);
    type_set_copy.emplace(ElementType::X);
    type_set_copy.insert(element.GetPseudo());
    element_initialized_cluster_hashmap_[element] = InitializeClusterHashMap(type_set_copy);
  }
  const auto num_atoms = reference_config.GetNumAtoms();
#pragma omp parallel for default(none)  shared(num_atoms, reference_config)
  {
    for (size_t i = 0; i < num_atoms; ++i) {
      for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
// #pragma omp critical
//         {
        site_bond_mapping_hashmap_[{i, j}] =
            GetClusterParametersMappingState(reference_config, {i, j});
        // }
      }
    }
  }
}

EnergyPredictor::~EnergyPredictor() = default;
std::pair<double, double> EnergyPredictor::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {

  return GetBarrierAndDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
std::pair<double, double> EnergyPredictor::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto &migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
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
    dE += base_theta_[i] * de_encode[i];
    e0 += theta[i] * e0_encode[i];
  }
  return {e0 + dE / 2, dE};
}
} // namespace pred
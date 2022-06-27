#include "StateChangePredictor.h"
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {
StateChangePredictor::StateChangePredictor(const std::string &predictor_filename,
                                           const cfg::Config &reference_config,
                                           std::set<Element> type_set)
    : type_set_(std::move(type_set)) {
  auto type_set_copy(type_set_);
  type_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(type_set_copy);

  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }
#pragma omp parallel for default(none) shared(reference_config, std::cout)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::cout << i << ' ' << j << std::endl;
      auto bond_mapping = GetClusterParametersMappingStateOfBond(reference_config, {i, j});
#pragma omp critical
      {
        site_bond_neighbors_hashmap_[{j, i}] = std::move(bond_mapping);
      }
    }
  }

}
} // namespace pred
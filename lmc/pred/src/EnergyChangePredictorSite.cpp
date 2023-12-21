#include "EnergyChangePredictorSite.h"
#include <omp.h>
#include <nlohmann/json.hpp>

namespace pred {
EnergyChangePredictorSite::EnergyChangePredictorSite(const std::string &predictor_filename,
                                                     const cfg::Config &reference_config,
                                                     std::set<Element> element_set)
    : element_set_(std::move(element_set)) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(element_set_copy);

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  nlohmann::json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      auto base_theta_json = parameters.at("theta");
      base_theta_ = {};
      for (const auto &theta : base_theta_json) {
        base_theta_.emplace_back(theta.get<double>());
      }
      // base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }
#pragma omp parallel for default(none) shared(reference_config)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    auto site_mapping = GetClusterParametersMappingStateSiteOf(reference_config, i);
#pragma omp critical
    {
      site_neighbors_hashmap_[i] = std::move(site_mapping);
    }
  }
}
EnergyChangePredictorSite::~EnergyChangePredictorSite() = default;
double EnergyChangePredictorSite::GetDeFromAtomIdSite(
    const cfg::Config &config, size_t atom_id, Element new_element) const {
  return GetDeFromLatticeIdSite(config, config.GetLatticeIdFromAtomId(atom_id), new_element);
}
double EnergyChangePredictorSite::GetDeFromLatticeIdSite(
    const cfg::Config &config, size_t lattice_id, Element new_element) const {
  const auto old_element = config.GetElementAtLatticeId(lattice_id);
  if (old_element == new_element) {
    return 0.0;
  }
  const auto mapping = site_neighbors_hashmap_.at(lattice_id);
  auto start_hashmap(initialized_cluster_hashmap_);
  auto end_hashmap(initialized_cluster_hashmap_);

  int label = 0;
  for (const auto &cluster_vector: mapping) {
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto lattice_id_in_cluster: cluster) {
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
        if (lattice_id_in_cluster == lattice_id) {
          element_vector_end.push_back(new_element);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
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
  static const std::vector<double>
      cluster_counter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto count_bond = static_cast<double>(end_hashmap.at(cluster))
        - static_cast<double>(start_hashmap.at(cluster));
    auto total_bond = cluster_counter[static_cast<size_t>(cluster.GetLabel())];
    de_encode.push_back(count_bond / total_bond);
  }
  double dE = 0;
  const size_t cluster_size = base_theta_.size();
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += base_theta_[i] * de_encode[i];
  }
  return dE;
}
} // pred
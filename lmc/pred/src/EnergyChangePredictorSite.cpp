#include "EnergyChangePredictorSite.h"
#include <omp.h>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <stdexcept>

namespace pred {
namespace {
const std::vector<double> kClusterCounter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
}
EnergyChangePredictorSite::EnergyChangePredictorSite(const std::string &predictor_filename,
                                                     const cfg::Config &reference_config,
                                                     std::set<Element> element_set)
    : element_set_(std::move(element_set)) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  const auto init_cluster_hashmap = InitializeClusterHashMap(element_set_copy);
  const std::map<cfg::ElementCluster, int> ordered(init_cluster_hashmap.begin(), init_cluster_hashmap.end());
  std::vector<double> cluster_total_bonds;
  cluster_total_bonds.reserve(ordered.size());
  for (const auto &entry : ordered) {
    cluster_total_bonds.push_back(kClusterCounter.at(static_cast<size_t>(entry.first.GetLabel())));
  }
  cluster_indexer_ = ClusterIndexer(ordered, std::move(cluster_total_bonds));

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  if (!ifs) {
    throw std::runtime_error("Cannot open " + predictor_filename);
  }
  nlohmann::json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = JsonToEigenVector(parameters.at("theta"));
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
  auto &start_counts = GetThreadLocalStartCountsBuffer();
  auto &end_counts = GetThreadLocalEndCountsBuffer();
  start_counts.assign(cluster_indexer_.Size(), 0);
  end_counts.assign(cluster_indexer_.Size(), 0);

  auto &element_vector_start = GetThreadLocalElementStartBuffer();
  auto &element_vector_end = GetThreadLocalElementEndBuffer();

  int label = 0;
  for (const auto &cluster_vector: mapping) {
    for (const auto &cluster: cluster_vector) {
      element_vector_start.clear();
      element_vector_end.clear();
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
      start_counts[cluster_indexer_.GetIndex(cfg::ElementCluster(label, element_vector_start))]++;
      end_counts[cluster_indexer_.GetIndex(cfg::ElementCluster(label, element_vector_end))]++;
    }
    label++;
  }

  auto &de_encode = GetThreadLocalDeEncodeBuffer();
  de_encode.resize(cluster_indexer_.Size());
  const auto &total_bonds = cluster_indexer_.GetTotalBonds();
  for (size_t idx = 0; idx < cluster_indexer_.Size(); ++idx) {
    de_encode[idx] = (static_cast<double>(end_counts[idx]) - static_cast<double>(start_counts[idx]))
        / total_bonds[idx];
  }
  const Eigen::Map<const Eigen::VectorXd> encode_vec(de_encode.data(), static_cast<Eigen::Index>(de_encode.size()));
  return base_theta_.dot(encode_vec);
}
} // pred

#include "EnergyChangePredictorPairSite.h"

#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <algorithm>
#include <stdexcept>

namespace pred {
namespace {
const std::vector<double> kClusterCounter{256, 1536, 768, 3072, 2048, 3072, 6144, 6144, 6144, 6144, 2048};
}
EnergyChangePredictorPairSite::EnergyChangePredictorPairSite(const std::string &predictor_filename,
                                                             const cfg::Config &reference_config,
                                                             std::set<Element> element_set)
    : element_set_(std::move(element_set)),
      site_mapping_state_(GetClusterParametersMappingStateSite(reference_config)) {
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
      auto base_theta_json = parameters.at("theta");
      base_theta_ = {};
      for (const auto &theta: base_theta_json) {
        base_theta_.emplace_back(theta.get<double>());
      }
      // base_theta_ = std::vector<double>(parameters.at("theta"));
    }
  }
#pragma omp parallel for default(none) shared(reference_config)
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    auto neighboring_hashset = reference_config.GetNeighborsLatticeIdSetOfSite(i);
    auto sorted_lattice_vector = GetSortedLatticeVectorStateOfSite(reference_config, i);
    std::vector<size_t> lattice_id_vector_state;
    std::transform(sorted_lattice_vector.begin(),
                   sorted_lattice_vector.end(),
                   std::back_inserter(lattice_id_vector_state),
                   [](const auto &lattice) {
                     return lattice.GetId();
                   });
#pragma omp critical
    {
      site_state_hashmap_[i] = std::move(lattice_id_vector_state);
      neighboring_sites_hashmap_[i] = std::move(neighboring_hashset);
    }
  }
}

EnergyChangePredictorPairSite::~EnergyChangePredictorPairSite() = default;

double EnergyChangePredictorPairSite::GetDeFromAtomIdPair(const cfg::Config &config,
                                                          const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  return GetDeFromLatticeIdPair(config,
                                {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
                                 config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}

double
EnergyChangePredictorPairSite::GetDeFromLatticeIdPair(const cfg::Config &config,
                                                      const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  if (config.GetElementAtLatticeId(lattice_id_jump_pair.first) ==
      config.GetElementAtLatticeId(lattice_id_jump_pair.second)) {
    return 0.0;
  }
  const auto &neighbors_of_lhs = neighboring_sites_hashmap_.at(lattice_id_jump_pair.first);
  if (neighbors_of_lhs.find(lattice_id_jump_pair.second) == neighbors_of_lhs.end()) {
    return GetDeFromLatticeIdPairWithoutCoupling(config, lattice_id_jump_pair);
  }
  return GetDeFromLatticeIdPairWithCoupling(config, lattice_id_jump_pair);
}

double EnergyChangePredictorPairSite::GetDeHelper(
    std::vector<int> &start_counts,
    std::vector<int> &end_counts) const {
  auto &de_encode = GetThreadLocalDoubleBuffer();
  de_encode.resize(cluster_indexer_.Size());
  const auto &total_bonds = cluster_indexer_.GetTotalBonds();
  for (size_t idx = 0; idx < cluster_indexer_.Size(); ++idx) {
    de_encode[idx] = (static_cast<double>(end_counts[idx]) - static_cast<double>(start_counts[idx]))
        / total_bonds[idx];
  }
  const Eigen::Map<const Eigen::VectorXd> theta_vec(base_theta_.data(), static_cast<Eigen::Index>(base_theta_.size()));
  const Eigen::Map<const Eigen::VectorXd> encode_vec(de_encode.data(), static_cast<Eigen::Index>(de_encode.size()));
  return theta_vec.dot(encode_vec);
}

double EnergyChangePredictorPairSite::GetDeFromLatticeIdPairWithCoupling(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
  const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  if (element_first == element_second) {
    return 0.0;
  }
  const auto mapping = GetClusterParametersMappingStatePairOf(config, lattice_id_jump_pair);
  auto &start_counts = GetThreadLocalIntPrimaryBuffer();
  auto &end_counts = GetThreadLocalIntSecondaryBuffer();
  start_counts.assign(cluster_indexer_.Size(), 0);
  end_counts.assign(cluster_indexer_.Size(), 0);
  auto &element_vector_start = GetThreadLocalElementPrimaryBuffer();
  auto &element_vector_end = GetThreadLocalElementSecondaryBuffer();
  int label = 0;
  for (const auto &cluster_vector: mapping) {
    for (const auto &cluster: cluster_vector) {
      element_vector_start.clear();
      element_vector_end.clear();
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto lattice_id: cluster) {
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id));
        if (lattice_id == lattice_id_jump_pair.first) {
          element_vector_end.push_back(element_second);
          continue;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          element_vector_end.push_back(element_first);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id));
      }
      start_counts[cluster_indexer_.GetIndex(cfg::ElementCluster(label, element_vector_start))]++;
      end_counts[cluster_indexer_.GetIndex(cfg::ElementCluster(label, element_vector_end))]++;
    }
    label++;
  }
  return GetDeHelper(start_counts, end_counts);
}

double EnergyChangePredictorPairSite::GetDeFromLatticeIdPairWithoutCoupling(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
  const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  if (element_first == element_second) {
    return 0.0;
  }
  auto dE1 = GetDeFromLatticeIdSite(config, lattice_id_jump_pair.first, element_second);
  auto dE2 = GetDeFromLatticeIdSite(config, lattice_id_jump_pair.second, element_first);
  return dE1 + dE2;
}

double EnergyChangePredictorPairSite::GetDeFromAtomIdSite(const cfg::Config &config,
                                                          size_t atom_id,
                                                          Element new_element) const {
  return GetDeFromLatticeIdSite(config, config.GetLatticeIdFromAtomId(atom_id), new_element);
}

double EnergyChangePredictorPairSite::GetDeFromLatticeIdSite(const cfg::Config &config,
                                                             size_t lattice_id,
                                                             Element new_element) const {
  const auto old_element = config.GetElementAtLatticeId(lattice_id);
  if (old_element == new_element) {
    return 0.0;
  }

  auto &start_counts = GetThreadLocalIntPrimaryBuffer();
  auto &end_counts = GetThreadLocalIntSecondaryBuffer();
  start_counts.assign(cluster_indexer_.Size(), 0);
  end_counts.assign(cluster_indexer_.Size(), 0);
  const auto &lattice_id_vector = site_state_hashmap_.at(lattice_id);
  auto &element_vector_start = GetThreadLocalElementPrimaryBuffer();
  auto &element_vector_end = GetThreadLocalElementSecondaryBuffer();

  int label = 0;
  for (const auto &cluster_vector: site_mapping_state_) {
    for (const auto &cluster: cluster_vector) {
      element_vector_start.clear();
      element_vector_end.clear();
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      for (auto index: cluster) {
        size_t lattice_id_in_cluster = lattice_id_vector[index];
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
        if (lattice_id_in_cluster == lattice_id) {
          element_vector_end.emplace_back(new_element);
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id_in_cluster));
      }
      start_counts[cluster_indexer_.GetIndex(cfg::ElementCluster(static_cast<int>(label), element_vector_start))]++;
      end_counts[cluster_indexer_.GetIndex(cfg::ElementCluster(static_cast<int>(label), element_vector_end))]++;
    }
    label++;
  }
  return GetDeHelper(start_counts, end_counts);
}
}    // namespace pred

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
      site_mapping_state_(FlattenClusterMapping(GetClusterParametersMappingStateSite(reference_config))) {
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
  CheckModelSize(predictor_filename, "Base.theta", static_cast<size_t>(base_theta_.size()), cluster_indexer_.Size());
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

double EnergyChangePredictorPairSite::GetDeHelper(const std::vector<int> &count_changes) const {
  auto &de_encode = GetThreadLocalDeEncodeBuffer();
  de_encode.resize(cluster_indexer_.Size());
  const auto &total_bonds = cluster_indexer_.GetTotalBonds();
  for (size_t idx = 0; idx < cluster_indexer_.Size(); ++idx) {
    de_encode[idx] = static_cast<double>(count_changes[idx]) / total_bonds[idx];
  }
  const Eigen::Map<const Eigen::VectorXd> encode_vec(de_encode.data(), static_cast<Eigen::Index>(de_encode.size()));
  return base_theta_.dot(encode_vec);
}

double EnergyChangePredictorPairSite::GetDeFromLatticeIdPairWithCoupling(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto element_first = config.GetElementAtLatticeId(lattice_id_jump_pair.first);
  const auto element_second = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  if (element_first == element_second) {
    return 0.0;
  }
  const auto mapping = GetClusterParametersMappingStatePairOf(config, lattice_id_jump_pair);
  // Number of clusters of every kind after the swap minus their number before it. The counts are small integers, so
  // the difference is exact.
  auto &count_changes = GetThreadLocalStartCountsBuffer();
  count_changes.assign(cluster_indexer_.Size(), 0);
  const auto code_first = cluster_indexer_.GetElementCode(element_first);
  const auto code_second = cluster_indexer_.GetElementCode(element_second);
  std::array<uint8_t, FlatClusterMapping::kMaxClusterSize> cluster_start{};
  std::array<uint8_t, FlatClusterMapping::kMaxClusterSize> cluster_end{};
  size_t label = 0;
  for (const auto &cluster_vector: mapping) {
    for (const auto &cluster: cluster_vector) {
      for (size_t i = 0; i < cluster.size(); ++i) {
        const auto lattice_id = cluster[i];
        const auto code = cluster_indexer_.GetElementCode(config.GetElementAtLatticeId(lattice_id));
        cluster_start[i] = code;
        if (lattice_id == lattice_id_jump_pair.first) {
          cluster_end[i] = code_second;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          cluster_end[i] = code_first;
        } else {
          cluster_end[i] = code;
        }
      }
      --count_changes[cluster_indexer_.GetIndex(label, cluster_start.data(), cluster.size())];
      ++count_changes[cluster_indexer_.GetIndex(label, cluster_end.data(), cluster.size())];
    }
    label++;
  }
  return GetDeHelper(count_changes);
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

  // Element codes of the neighbourhood before and after the site has changed its element.
  const auto &lattice_id_vector = site_state_hashmap_.at(lattice_id);
  const auto new_code = cluster_indexer_.GetElementCode(new_element);
  auto &start_codes = GetThreadLocalStartCodesBuffer();
  auto &end_codes = GetThreadLocalEndCodesBuffer();
  start_codes.resize(lattice_id_vector.size());
  end_codes.resize(lattice_id_vector.size());
  for (size_t i = 0; i < lattice_id_vector.size(); ++i) {
    const size_t lattice_id_in_cluster = lattice_id_vector[i];
    const auto code = cluster_indexer_.GetElementCode(config.GetElementAtLatticeId(lattice_id_in_cluster));
    start_codes[i] = code;
    end_codes[i] = lattice_id_in_cluster == lattice_id ? new_code : code;
  }

  auto &count_changes = GetThreadLocalStartCountsBuffer();
  count_changes.assign(cluster_indexer_.Size(), 0);
  AddClusterCountChanges(site_mapping_state_, cluster_indexer_, start_codes, end_codes, count_changes);
  return GetDeHelper(count_changes);
}
}    // namespace pred

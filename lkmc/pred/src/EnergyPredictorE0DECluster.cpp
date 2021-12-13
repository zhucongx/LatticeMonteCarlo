#include "EnergyPredictorE0DECluster.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {

static std::unordered_map<ElementCluster,
                          size_t,
                          boost::hash<ElementCluster> > InitializeClusterHashMap(
    const std::set<Element> &type_set) {
  std::unordered_map<ElementCluster, size_t,
                     boost::hash<ElementCluster> > initialized_cluster_hashmap;
  for (const auto &element1: type_set) {
    initialized_cluster_hashmap[ElementCluster(0, element1)] = 0;
  }
  for (size_t label = 1; label <= 3; ++label) {
    for (const auto &element1: type_set) {
      for (const auto &element2: type_set) {
        initialized_cluster_hashmap[ElementCluster(label, element1, element2)] = 0;
      }
    }
  }
  for (size_t label = 4; label < 11; ++label) {
    for (const auto &element1: type_set) {
      for (const auto &element2: type_set) {
        for (const auto &element3: type_set) {
          initialized_cluster_hashmap[ElementCluster(label, element1, element2, element3)] = 0;
        }
      }
    }
  }
  return initialized_cluster_hashmap;
}

static std::vector<double> GetClusterChange(
    const std::vector<Element> &encode,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping_periodic,
    std::unordered_map<ElementCluster, size_t,
                       boost::hash<ElementCluster> > initialized_cluster_hashmap) {
  for (size_t label = 0; label < 11; ++label) {
    for (const auto &cluster: cluster_mapping_periodic[label]) {
      std::vector<Element> cluster_element_vector;
      cluster_element_vector.reserve(cluster.size());
      for (auto index: cluster) {
        cluster_element_vector.push_back(encode[index]);
      }
      initialized_cluster_hashmap[ElementCluster{label, cluster_element_vector}]++;
    }
  }

  for (auto &[element_cluster, count]: initialized_cluster_hashmap) {
    if (element_cluster.GetClusterSize() == 2) { count /= 2; }
    else if (element_cluster.GetClusterSize() == 3) { count /= 6; }
  }

  std::map<ElementCluster, int>
      ordered(initialized_cluster_hashmap.begin(), initialized_cluster_hashmap.end());
  std::vector<double> res;
  res.reserve(ordered.size());
  for (const auto &bond_count: ordered) {
    res.push_back(bond_count.second);
  }
  return res;
}

EnergyPredictorE0DECluster::EnergyPredictorE0DECluster(const std::string &predictor_filename,
                                                       const cfg::Config &reference_config,
                                                       const std::set<Element> &type_set)
    : EnergyPredictor(type_set),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_periodic_(GetAverageClusterParametersMappingPeriodic(reference_config)),
      initialized_cluster_hashmap_(InitializeClusterHashMap(type_set_)) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Cluster") {
      theta_cluster_ = std::vector<double>(parameters.at("theta"));
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

EnergyPredictorE0DECluster::~EnergyPredictorE0DECluster() = default;

double EnergyPredictorE0DECluster::GetE0(const cfg::Config &config,
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
double EnergyPredictorE0DECluster::GetDE(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    Element migration_element) const {
  auto lattice_id_vector_periodic_start =
      site_cluster_periodic_hashmap_.at(lattice_id_jump_pair.first);

  std::vector<Element> element_vector_start{}, element_vector_end{};
  element_vector_start.reserve(lattice_id_vector_periodic_start.size());
  for (auto index: lattice_id_vector_periodic_start) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementType::X) {
      element_vector_start.push_back(migration_element);
      continue;
    }
    element_vector_start.push_back(this_element);
  }
  auto start_encode = GetClusterChange(element_vector_start,
                                       mapping_periodic_,
                                       initialized_cluster_hashmap_);

  auto lattice_id_vector_periodic_end =
      site_cluster_periodic_hashmap_.at(lattice_id_jump_pair.second);
  element_vector_end.reserve(lattice_id_vector_periodic_end.size());
  for (auto index: lattice_id_vector_periodic_end) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementType::X) {
      element_vector_end.push_back(migration_element);
      continue;
    }
    element_vector_end.push_back(this_element);
  }
  auto end_encode = GetClusterChange(element_vector_end,
                                     mapping_periodic_,
                                     initialized_cluster_hashmap_);

  // for (size_t i = 0; i < end_encode.size(); ++i) {
  //   std::cout << end_encode[i] - start_encode[i] << ", ";
  // }
  // std::cout << "\n";
  const size_t bond_size = theta_cluster_.size();
  double dE = 0;
  for (size_t i = 0; i < bond_size; ++i) {
    dE += theta_cluster_[i] * (end_encode[i] - start_encode[i]);
  }
  return dE;
}
std::pair<double, double> EnergyPredictorE0DECluster::GetBarrierAndDiffFromLatticeIdPair(
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

  auto dE = GetDE(config, lattice_id_jump_pair, migration_element);

  return {e0 + dE / 2, dE};
}

} // namespace pred
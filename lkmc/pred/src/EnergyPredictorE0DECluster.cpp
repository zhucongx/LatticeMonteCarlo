#include "EnergyPredictorE0DECluster.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace pred {


static std::unordered_map<ElementCluster, int,
                          boost::hash<ElementCluster> > GetClusterEncode(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    std::unordered_map<ElementCluster, int,
                       boost::hash<ElementCluster> > initialized_cluster_hashmap) {
  auto lattice_id_set = GetNeighborsLatticeIdSetOfJumpPair(config, lattice_id_jump_pair);
  for (auto lattice1_id: lattice_id_set) {
    const Element type1 = config.GetElementAtLatticeId(lattice1_id);
    // 0 singlets not changing
    // plus new pairs
    for (auto lattice2_id: config.GetFirstNeighborsAdjacencyList()[lattice1_id]) {
      if (lattice_id_set.find(lattice2_id) == lattice_id_set.end()) { continue; }
      // 1 first pair
      Element type2 = config.GetElementAtLatticeId(lattice2_id);
      initialized_cluster_hashmap[ElementCluster{1, type1, type2}]++;
      for (auto lattice3_id: config.GetFirstNeighborsAdjacencyList()[lattice2_id]) {
        if (lattice_id_set.find(lattice3_id) == lattice_id_set.end()) { continue; }
        Element type3 = config.GetElementAtLatticeId(lattice3_id);
        if (std::find(config.GetFirstNeighborsAdjacencyList()[lattice1_id].begin(),
                      config.GetFirstNeighborsAdjacencyList()[lattice1_id].end(),
                      lattice3_id)
            != config.GetFirstNeighborsAdjacencyList()[lattice1_id].end()) {
          // 4 first-first-first triplet
          initialized_cluster_hashmap[ElementCluster{4, type1, type2, type3}]++;
        } else if (std::find(config.GetSecondNeighborsAdjacencyList()[lattice1_id].begin(),
                             config.GetSecondNeighborsAdjacencyList()[lattice1_id].end(),
                             lattice3_id)
            != config.GetSecondNeighborsAdjacencyList()[lattice1_id].end()) {
          // 5 first-first-second triplet
          initialized_cluster_hashmap[ElementCluster{5, type1, type2, type3}]++;

        } else if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_id].begin(),
                             config.GetThirdNeighborsAdjacencyList()[lattice1_id].end(),
                             lattice3_id)
            != config.GetThirdNeighborsAdjacencyList()[lattice1_id].end()) {
          // 6 first-first-third triplet
          initialized_cluster_hashmap[ElementCluster{6, type1, type2, type3}]++;
        }
      }
      for (auto lattice3_id: config.GetSecondNeighborsAdjacencyList()[lattice2_id]) {
        if (lattice_id_set.find(lattice3_id) == lattice_id_set.end()) { continue; }
        Element type3 = config.GetElementAtLatticeId(lattice3_id);
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_id].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_id].end(),
                      lattice3_id)
            != config.GetThirdNeighborsAdjacencyList()[lattice1_id].end()) {
          // 7 first-second-third triplet
          initialized_cluster_hashmap[ElementCluster{7, type1, type2, type3}]++;
        }
      }
      for (auto lattice3_id: config.GetThirdNeighborsAdjacencyList()[lattice2_id]) {
        if (lattice_id_set.find(lattice3_id) == lattice_id_set.end()) { continue; }
        Element type3 = config.GetElementAtLatticeId(lattice3_id);
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_id].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_id].end(),
                      lattice3_id)
            != config.GetThirdNeighborsAdjacencyList()[lattice1_id].end()) {
          // 8 first-third-third triplet
          initialized_cluster_hashmap[ElementCluster{8, type1, type2, type3}]++;
        }
      }
    }
    for (auto lattice2_id: config.GetSecondNeighborsAdjacencyList()[lattice1_id]) {
      if (lattice_id_set.find(lattice2_id) == lattice_id_set.end()) { continue; }
      // 2 second pair
      Element type2 = config.GetElementAtLatticeId(lattice2_id);
      initialized_cluster_hashmap[ElementCluster{2, type1, type2}]++;
      for (auto lattice3_id: config.GetThirdNeighborsAdjacencyList()[lattice2_id]) {
        if (lattice_id_set.find(lattice3_id) == lattice_id_set.end()) { continue; }
        Element type3 = config.GetElementAtLatticeId(lattice3_id);
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_id].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_id].end(),
                      lattice3_id)
            != config.GetThirdNeighborsAdjacencyList()[lattice1_id].end()) {
          // 9 second-third-third triplet
          initialized_cluster_hashmap[ElementCluster{9, type1, type2, type3}]++;
        }
      }
    }
    for (auto lattice2_id: config.GetThirdNeighborsAdjacencyList()[lattice1_id]) {
      if (lattice_id_set.find(lattice2_id) == lattice_id_set.end()) { continue; }
      // 3 second pair
      Element type2 = config.GetElementAtLatticeId(lattice2_id);
      initialized_cluster_hashmap[ElementCluster{3, type1, type2}]++;
      for (auto lattice3_id: config.GetThirdNeighborsAdjacencyList()[lattice2_id]) {
        if (lattice_id_set.find(lattice3_id) == lattice_id_set.end()) { continue; }
        Element type3 = config.GetElementAtLatticeId(lattice3_id);
        if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice1_id].begin(),
                      config.GetThirdNeighborsAdjacencyList()[lattice1_id].end(),
                      lattice3_id)
            != config.GetThirdNeighborsAdjacencyList()[lattice1_id].end()) {
          // 10 third-third-third triplet
          initialized_cluster_hashmap[ElementCluster{10, type1, type2, type3}]++;
        }
      }
    }
  }
  return initialized_cluster_hashmap;
}
static std::vector<double> GetClusterChange(
    cfg::Config config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    std::unordered_map<ElementCluster, int,
                       boost::hash<ElementCluster> > initialized_cluster_hashmap) {

  const auto start_cluster_encode =
      GetClusterEncode(config, lattice_id_jump_pair, initialized_cluster_hashmap);
  config.LatticeJump(lattice_id_jump_pair);
  const auto end_cluster_encode =
      GetClusterEncode(config, lattice_id_jump_pair, initialized_cluster_hashmap);

  std::map<ElementCluster, int>
      ordered(initialized_cluster_hashmap.begin(), initialized_cluster_hashmap.end());
  std::vector<double> res;
  res.reserve(ordered.size());
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto change =
        static_cast<double>(end_cluster_encode.at(cluster) - start_cluster_encode.at(cluster));
    double total_bond;
    switch (cluster.GetLabel()) {
      case 0:total_bond = 256;
        break;
      case 1: total_bond = 3072;
        break;
      case 2: total_bond = 1536;
        break;
      case 3: total_bond = 6144;
        break;
      case 4: total_bond = 12288;
        break;
      case 5: total_bond = 6144;
        break;
      case 6: total_bond = 12288;
        break;
      case 7: total_bond = 6144;
        break;
      case 8: total_bond = 12288;
        break;
      case 9: total_bond = 12288;
        break;
      case 10: total_bond = 12288;
        break;
    }
    res.push_back(change / total_bond);
  }
  return res;
}

EnergyPredictorE0DECluster::EnergyPredictorE0DECluster(const std::string &predictor_filename,
                                                       const cfg::Config &reference_config,
                                                       const std::set<Element> &type_set)
    : EnergyPredictor(type_set),
      one_hot_encode_hash_map_(GetOneHotEncodeHashmap(type_set)),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)) {
  auto type_set_copy = type_set;
  type_set_copy.emplace("X");
  type_set_copy.emplace("pAl");
  initialized_cluster_hashmap_ = InitializeClusterHashMap(type_set_copy);
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Cluster") {
      theta_cluster_ = std::vector<double>(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[Element(element)] = ParametersE0{
        parameters.at("e0_mu_x"),
        parameters.at("e0_sigma_x"),
        parameters.at("e0_mu_y"),
        parameters.at("e0_sigma_y"),
        parameters.at("e0_U"),
        parameters.at("e0_theta"),
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
double EnergyPredictorE0DECluster::GetDE(const cfg::Config &config,
                                         const std::pair<size_t,
                                                         size_t> &lattice_id_jump_pair) const {
  auto cluster_change_vector =
      GetClusterChange(config, lattice_id_jump_pair, initialized_cluster_hashmap_);
  const size_t cluster_size = theta_cluster_.size();
  double dE = 0;
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += theta_cluster_[i] * cluster_change_vector[i];
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

  auto dE = GetDE(config, lattice_id_jump_pair);
  return {e0 + dE / 2, dE};
}

} // namespace pred
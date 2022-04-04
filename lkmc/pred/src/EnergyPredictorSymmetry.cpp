#include "EnergyPredictorSymmetry.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace pred {
using Singlet_MMM_t = cfg::LatticeClusterMMM<1>;
using Pair_MMM_t = cfg::LatticeClusterMMM<2>;
using Singlet_MM2_t = cfg::LatticeClusterMM2<1>;
using Pair_MM2_t = cfg::LatticeClusterMM2<2>;

static std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<Element> &type_set) {
  size_t type_size = type_set.size();
  std::unordered_map<std::string, std::vector<double> > encode_dict;

  size_t ct1 = 0;
  for (const auto &element: type_set) {
    std::vector<double> element_encode(type_size, 0);
    element_encode[ct1] = 1.0;
    encode_dict[element.GetString()] = element_encode;
    ++ct1;
  }

  size_t num_pairs = type_size * type_size;
  size_t ct2 = 0;
  for (const auto &element1: type_set) {
    for (const auto &element2: type_set) {
      std::vector<double> element_encode(num_pairs, 0);
      element_encode[ct2] = 1.0;
      encode_dict[element1.GetString() + element2.GetString()] = element_encode;
      ++ct2;
    }
  }

  size_t num_triplets = type_size * type_size * type_size;
  size_t ct3 = 0;
  for (const auto &element1: type_set) {
    for (const auto &element2: type_set) {
      for (const auto &element3: type_set) {
        std::vector<double> element_encode(num_triplets, 0);
        element_encode[ct3] = 1.0;
        encode_dict[element1.GetString() + element2.GetString() + element3.GetString()] =
            element_encode;
        ++ct3;
      }
    }
  }
  return encode_dict;
}

static bool LatticeSortCompareMMM(const cfg::Lattice &lhs,
                                  const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double diff_norm = Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
  if (diff_norm < -kEpsilon)
    return true;
  if (diff_norm > kEpsilon)
    return false;
  const double diff_x_sym = std::abs(relative_position_lhs[kXDimension] - 0.5)
      - std::abs(relative_position_rhs[kXDimension] - 0.5);
  return diff_x_sym < -kEpsilon;
}
static bool LatticeSortCompareMM2(const cfg::Lattice &lhs,
                                  const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double diff_norm = Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
  if (diff_norm < -kEpsilon)
    return true;
  if (diff_norm > kEpsilon)
    return false;
  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  return diff_x < -kEpsilon;
}

template<size_t DataSize>
static bool IsClusterSmallerSymmetricallyMMM(const cfg::LatticeClusterMMM<DataSize> &lhs,
                                             const cfg::LatticeClusterMMM<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    const auto &lhs_lattice = lhs.GetLatticeAt(i);
    const auto &rhs_lattice = rhs.GetLatticeAt(i);
    if (LatticeSortCompareMMM(lhs_lattice, rhs_lattice)) { return true; }
    if (LatticeSortCompareMMM(rhs_lattice, lhs_lattice)) { return false; }
  }
  // if it reaches here, it means that the clusters are same symmetrically. Returns false.
  return false;
}
template<size_t DataSize>
static bool IsClusterSmallerSymmetricallyMM2(const cfg::LatticeClusterMM2<DataSize> &lhs,
                                             const cfg::LatticeClusterMM2<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    const auto &lhs_lattice = lhs.GetLatticeAt(i);
    const auto &rhs_lattice = rhs.GetLatticeAt(i);
    if (LatticeSortCompareMM2(lhs_lattice, rhs_lattice))
      return true;
    if (LatticeSortCompareMM2(rhs_lattice, lhs_lattice))
      return false;
  }
  // if it reaches here, it means that the clusters are same symmetrically. Returns false.
  return false;
}

// Returns forward and backward sorted lattice lists
static std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMMM(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  // The number of first-, second-, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = 60;
  auto lattice_id_hashset =
      GetNeighborsLatticeIdSetOfJumpPair(config, lattice_id_jump_pair);
  const auto move_distance = Vector_t{0.5, 0.5, 0.5}
      - GetLatticePairCenter(config, lattice_id_jump_pair);
  std::vector<cfg::Lattice> lattice_list;
  lattice_list.reserve(kNumOfSites);
  for (const auto id: lattice_id_hashset) {
    cfg::Lattice lattice = config.GetLatticeVector()[id];
    // move to center
    auto relative_position = lattice.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);
    lattice.SetRelativePosition(relative_position);
    if (lattice.GetId() == lattice_id_jump_pair.first
        || lattice.GetId() == lattice_id_jump_pair.second)
      continue;
    lattice_list.push_back(lattice);
  }
  RotateLatticeVector(lattice_list,
                      GetLatticePairRotationMatrix(config, lattice_id_jump_pair));
  //sort using mmm group point (mm2 if x mirror not applied)
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              if (LatticeSortCompareMMM(lhs, rhs)) return true;
              if (LatticeSortCompareMMM(rhs, lhs)) return false;
              const auto &relative_position_lhs = lhs.GetRelativePosition();
              const auto &relative_position_rhs = rhs.GetRelativePosition();
              const double diff_x =
                  relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
              if (diff_x < -kEpsilon) { return true; }
              if (diff_x > kEpsilon) { return false; }
              const double diff_y =
                  relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
              if (diff_y < -kEpsilon) { return true; }
              if (diff_y > kEpsilon) { return false; }
              const double diff_z =
                  relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
              if (diff_z < -kEpsilon) { return true; }
              if (diff_z > kEpsilon) { return false; }
              return lhs.GetId() < rhs.GetId();
            });
  return lattice_list;
}
static std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMM2(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  // The number of first-, second-, and third-nearest neighbors of the jump pairs
  constexpr size_t kNumOfSites = 60;
  auto lattice_id_hashset =
      GetNeighborsLatticeIdSetOfJumpPair(config, lattice_id_jump_pair);
  const auto move_distance = Vector_t{0.5, 0.5, 0.5}
      - GetLatticePairCenter(config, lattice_id_jump_pair);
  std::vector<cfg::Lattice> lattice_list;
  lattice_list.reserve(kNumOfSites);
  for (const auto id: lattice_id_hashset) {
    cfg::Lattice lattice = config.GetLatticeVector()[id];
    // move to center
    auto relative_position = lattice.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);
    lattice.SetRelativePosition(relative_position);
    if (lattice.GetId() == lattice_id_jump_pair.first
        || lattice.GetId() == lattice_id_jump_pair.second)
      continue;
    lattice_list.push_back(lattice);
  }
  RotateLatticeVector(lattice_list,
                      GetLatticePairRotationMatrix(config, lattice_id_jump_pair));
  //sort using mmm group point (mm2 if x mirror not applied)
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              if (LatticeSortCompareMM2(lhs, rhs)) return true;
              if (LatticeSortCompareMM2(rhs, lhs)) return false;
              const auto &relative_position_lhs = lhs.GetRelativePosition();
              const auto &relative_position_rhs = rhs.GetRelativePosition();
              const double diff_x =
                  relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
              if (diff_x < -kEpsilon) { return true; }
              if (diff_x > kEpsilon) { return false; }
              const double diff_y =
                  relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
              if (diff_y < -kEpsilon) { return true; }
              if (diff_y > kEpsilon) { return false; }
              const double diff_z =
                  relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
              if (diff_z < -kEpsilon) { return true; }
              if (diff_z > kEpsilon) { return false; }
              return lhs.GetId() < rhs.GetId();
            });
  return lattice_list;
}

template<size_t DataSize>
static void GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
    std::vector<cfg::LatticeClusterMMM<DataSize> > &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  // sort clusters
  std::sort(cluster_vector.begin(), cluster_vector.end(),
            [](const auto &lhs, const auto &rhs) {
              return IsClusterSmallerSymmetricallyMMM(lhs, rhs);
            });
  // start to point at Cluster in the first range
  typename std::vector<cfg::LatticeClusterMMM<DataSize> >::const_iterator lower_it, upper_it;
  lower_it = cluster_vector.cbegin();
  do {
    upper_it = std::upper_bound(lower_it, cluster_vector.cend(),
                                *lower_it,
                                [](const auto &lhs, const auto &rhs) {
                                  return IsClusterSmallerSymmetricallyMMM(lhs, rhs);
                                });
    std::vector<std::vector<size_t> > cluster_index_vector;
    for (auto it = lower_it; it != upper_it; ++it) {
      auto cluster_index = it->GetIndexVector();
      cluster_index_vector.push_back(cluster_index);
    }
    cluster_mapping.push_back(cluster_index_vector);
    // update to next range
    lower_it = upper_it;
  } while (upper_it != cluster_vector.cend());
}
template<size_t DataSize>
static void GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
    std::vector<cfg::LatticeClusterMM2<DataSize> > &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {
  // sort clusters
  std::sort(cluster_vector.begin(), cluster_vector.end(),
            [](const auto &lhs, const auto &rhs) {
              return IsClusterSmallerSymmetricallyMM2(lhs, rhs);
            });
  // start to point at Cluster in the first range
  typename std::vector<cfg::LatticeClusterMM2<DataSize> >::const_iterator lower_it, upper_it;
  lower_it = cluster_vector.cbegin();
  do {
    upper_it = std::upper_bound(lower_it, cluster_vector.cend(),
                                *lower_it,
                                [](const auto &lhs, const auto &rhs) {
                                  return IsClusterSmallerSymmetricallyMM2(lhs, rhs);
                                });
    std::vector<std::vector<size_t> > cluster_index_vector;
    for (auto it = lower_it; it != upper_it; ++it) {
      auto cluster_index = it->GetIndexVector();
      cluster_index_vector.push_back(cluster_index);
    }
    cluster_mapping.push_back(cluster_index_vector);
    // update to next range
    lower_it = upper_it;
  } while (upper_it != cluster_vector.cend());
}

static std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMM2(
    const cfg::Config &config) {
  const std::pair<size_t, size_t>
      lattice_id_jump_pair = {0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto lattice_vector = GetSymmetricallySortedLatticeVectorMM2(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  /// singlets
  std::vector<Singlet_MM2_t> singlet_vector;
  /// pairs
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > first_pair_set;
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > second_pair_set;
  std::unordered_set<Pair_MM2_t, boost::hash<Pair_MM2_t> > third_pair_set;
  for (size_t index1 = 0; index1 < lattice_vector.size(); ++index1) {
    cfg::Lattice lattice1(lattice_vector[index1]);
    const size_t id1 = lattice1.GetId();
    lattice1.SetId(index1);
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    for (size_t index2 = 0; index2 < index1; ++index2) {
      cfg::Lattice lattice2(lattice_vector[index2]);
      const size_t id2 = lattice2.GetId();
      lattice2.SetId(index2);
      switch (FindDistanceLabelBetweenLattice(id1, id2, config)) {
        case 1: first_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 2: second_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 3: third_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
      }
    }
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::move(singlet_vector), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::vector<Pair_MM2_t>(first_pair_set.begin(),
                              first_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::vector<Pair_MM2_t>(second_pair_set.begin(),
                              second_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMM2(
      std::vector<Pair_MM2_t>(third_pair_set.begin(),
                              third_pair_set.end()), cluster_mapping);

  return cluster_mapping;
}
static std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMMM(
    const cfg::Config &config) {
  const std::pair<size_t, size_t>
      lattice_id_jump_pair = {0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto lattice_vector = GetSymmetricallySortedLatticeVectorMMM(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  /// singlets
  std::vector<Singlet_MMM_t> singlet_vector;
  /// pairs
  std::unordered_set<Pair_MMM_t, boost::hash<Pair_MMM_t> > first_pair_set;
  std::unordered_set<Pair_MMM_t, boost::hash<Pair_MMM_t> > second_pair_set;
  std::unordered_set<Pair_MMM_t, boost::hash<Pair_MMM_t> > third_pair_set;
  for (size_t index1 = 0; index1 < lattice_vector.size(); ++index1) {
    cfg::Lattice lattice1(lattice_vector[index1]);
    const size_t id1 = lattice1.GetId();
    lattice1.SetId(index1);
    singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    for (size_t index2 = 0; index2 < index1; ++index2) {
      cfg::Lattice lattice2(lattice_vector[index2]);
      const size_t id2 = lattice2.GetId();
      lattice2.SetId(index2);
      switch (FindDistanceLabelBetweenLattice(id1, id2, config)) {
        case 1: first_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 2: second_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
        case 3: third_pair_set.emplace(std::array<cfg::Lattice, 2>{lattice1, lattice2});
          break;
      }
    }
  }
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::move(singlet_vector), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::vector<Pair_MMM_t>(first_pair_set.begin(),
                              first_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::vector<Pair_MMM_t>(second_pair_set.begin(),
                              second_pair_set.end()), cluster_mapping);
  GetAverageParametersMappingFromLatticeClusterVectorHelperMMM(
      std::vector<Pair_MMM_t>(third_pair_set.begin(),
                              third_pair_set.end()), cluster_mapping);
  return cluster_mapping;
}

static std::vector<double> GetOneHotParametersFromMap(
    const std::vector<Element> &encode,
    const std::unordered_map<std::string, std::vector<double> > &one_hot_encode_hashmap,
    size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping) {

  std::vector<double> res_encode;
  res_encode.reserve(354); // Todo check this magic number
  for (const auto &cluster_vector: cluster_mapping) {
    std::vector<double> sum_of_list(
        static_cast<size_t>(std::pow(num_of_elements, cluster_vector[0].size())), 0);
    for (const auto &cluster: cluster_vector) {
      std::string cluster_type;
      for (auto index: cluster) {
        cluster_type += encode[index].GetString();
      }
      const auto &cluster_one_hot_encode = one_hot_encode_hashmap.at(cluster_type);
      std::transform(sum_of_list.begin(), sum_of_list.end(),
                     cluster_one_hot_encode.begin(),
                     sum_of_list.begin(),
                     std::plus<>());
    }
    auto cluster_vector_size = static_cast<double>( cluster_vector.size());
    std::for_each(sum_of_list.begin(),
                  sum_of_list.end(),
                  [cluster_vector_size](auto &n) { n /= cluster_vector_size; });

    std::move(sum_of_list.begin(), sum_of_list.end(), std::back_inserter(res_encode));
  }
  return res_encode;
}

pred::EnergyPredictorSymmetry::EnergyPredictorSymmetry(const std::string &predictor_filename,
                                                       const cfg::Config &reference_config,
                                                       const std::set<Element> &type_set)
    : one_hot_encode_hash_map_(GetOneHotEncodeHashmap(type_set)),
      mapping_mmm_(GetAverageClusterParametersMappingMMM(reference_config)),
      mapping_mm2_(GetAverageClusterParametersMappingMM2(reference_config)),
      type_set_(type_set) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters]: all_parameters.items()) {
    e0_element_parameters_hashmap_[Element(element)] = ParametersE0{
        parameters.at("e0_mu_x"),
        parameters.at("e0_sigma_x"),
        parameters.at("e0_mu_y"),
        parameters.at("e0_sigma_y"),
        parameters.at("e0_U"),
        parameters.at("e0_theta"),
    };
    dE_element_parameters_hashmap_[Element(element)] = ParametersDE{
        parameters.at("dE_mu_x"),
        parameters.at("dE_sigma_x"),
        parameters.at("dE_mu_y"),
        parameters.at("dE_sigma_y"),
        parameters.at("dE_U"),
        parameters.at("dE_theta"),
    };
  }
  for (size_t i = 0; i < reference_config.GetNumAtoms(); ++i) {
    for (auto j: reference_config.GetFirstNeighborsAdjacencyList()[i]) {
      auto sorted_lattice_vector_mmm =
          GetSymmetricallySortedLatticeVectorMMM(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mmm;
      std::transform(sorted_lattice_vector_mmm.begin(), sorted_lattice_vector_mmm.end(),
                     std::back_inserter(lattice_id_vector_mmm),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_mmm_hashmap_[{i, j}] = lattice_id_vector_mmm;

      auto sorted_lattice_vector_mm2 =
          GetSymmetricallySortedLatticeVectorMM2(reference_config, {i, j});
      std::vector<size_t> lattice_id_vector_mm2;
      std::transform(sorted_lattice_vector_mm2.begin(), sorted_lattice_vector_mm2.end(),
                     std::back_inserter(lattice_id_vector_mm2),
                     [](const auto &lattice) { return lattice.GetId(); });
      site_bond_cluster_mm2_hashmap_[{i, j}] = lattice_id_vector_mm2;
    }
  }
}
EnergyPredictorSymmetry::~EnergyPredictorSymmetry() = default;
std::pair<double, double> EnergyPredictorSymmetry::GetBarrierAndDiffFromAtomIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {

  return GetBarrierAndDiffFromLatticeIdPair(
      config,
      {config.GetLatticeIdFromAtomId(atom_id_jump_pair.first),
       config.GetLatticeIdFromAtomId(atom_id_jump_pair.second)});
}
std::pair<double, double> EnergyPredictorSymmetry::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  auto e01 = GetE0(config,
                   lattice_id_jump_pair,
                   migration_element);
  auto e02 = GetE0(config,
                   {lattice_id_jump_pair.second, lattice_id_jump_pair.first},
                   migration_element);
  std::cerr << e01 << ", " << e02 << std::endl;
  auto e0 = (e01 + e02) / 2;
  auto ef = GetE(config,
                 lattice_id_jump_pair,
                 migration_element);
  auto eb = GetE(config,
                 {lattice_id_jump_pair.second, lattice_id_jump_pair.first},
                 migration_element);
  auto dE = eb - ef;
  return {e0 + dE / 2, dE};
}
double EnergyPredictorSymmetry::GetE0(const cfg::Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                      Element migration_element) const {
  auto lattice_id_vector_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mmm.size());
  for (auto index: lattice_id_vector_mmm) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementName::X) {
      ele_vector.push_back(migration_element);
      continue;
    }
    ele_vector.push_back(this_element);
  }
  auto encode_mmm = pred::GetOneHotParametersFromMap(ele_vector,
                                                     one_hot_encode_hash_map_,
                                                     type_set_.size(),
                                                     mapping_mmm_);
  const auto &element_parameters = e0_element_parameters_hashmap_.at(migration_element);

  const auto &mu_x = element_parameters.mu_x;
  const auto &sigma_x = element_parameters.sigma_x;
  const auto mu_y = element_parameters.mu_y;
  const auto sigma_y = element_parameters.sigma_y;

  const auto &U = element_parameters.U;
  const auto &theta = element_parameters.theta;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode_mmm.at(i) -= mu_x.at(i);
    encode_mmm.at(i) /= sigma_x.at(i);
  }
  double e0 = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mmm.at(i) * U.at(j).at(i);
    }
    e0 += pca_dot * theta.at(j);
  }
  e0 *= sigma_y;
  e0 += mu_y;
  return e0;
}
double EnergyPredictorSymmetry::GetE(const cfg::Config &config,
                                     const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                     Element migration_element) const {
  auto lattice_id_vector_mm2 = site_bond_cluster_mm2_hashmap_.at(lattice_id_jump_pair);
  std::vector<Element> ele_vector{};
  ele_vector.reserve(lattice_id_vector_mm2.size());
  for (auto index: lattice_id_vector_mm2) {
    auto this_element = config.GetElementAtLatticeId(index);
    if (this_element == ElementName::X) {
      ele_vector.push_back(migration_element);
      continue;
    }
    ele_vector.push_back(this_element);
  }
  auto encode_mm2 = pred::GetOneHotParametersFromMap(ele_vector,
                                                     one_hot_encode_hash_map_,
                                                     type_set_.size(),
                                                     mapping_mm2_);
  const auto &element_parameters = dE_element_parameters_hashmap_.at(migration_element);

  const auto &mu_x = element_parameters.mu_x;
  const auto &sigma_x = element_parameters.sigma_x;
  const auto mu_y = element_parameters.mu_y;
  const auto sigma_y = element_parameters.sigma_y;

  const auto &U = element_parameters.U;
  const auto &theta = element_parameters.theta;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    encode_mm2.at(i) -= mu_x.at(i);
    encode_mm2.at(i) /= sigma_x.at(i);
  }
  double e = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += encode_mm2.at(i) * U.at(j).at(i);
    }
    e += pca_dot * theta.at(j);
  }
  e *= sigma_y;
  e += mu_y;
  return e;
}
} // namespace pred

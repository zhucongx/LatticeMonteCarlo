#include "EnergyPredictor.h"
#include <utility>
#include <boost/range/combine.hpp>
#include <boost/functional/hash.hpp>

#include <omp.h>
#include <nlohmann/json.hpp>
#include "LatticeCluster.hpp"

using json = nlohmann::json;

namespace pred {
using Singlet_State_t = cfg::LatticeCluster<1>;
using Pair_State_t = cfg::LatticeCluster<2>;
using Triplet_State_t = cfg::LatticeCluster<3>;

static bool LatticeSortCompare(const cfg::Lattice &lhs,
                               const cfg::Lattice &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();
  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  const double diff_y = relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
  const double diff_z = relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
  if (diff_x < -kEpsilon) { return true; }
  if (diff_x > kEpsilon) { return false; }
  if (diff_y < -kEpsilon) { return true; }
  if (diff_y > kEpsilon) { return false; }
  return diff_z < -kEpsilon;
}

// Returns forward and backward sorted lattice lists
static std::vector<cfg::Lattice> GetSortedLatticeVectorState(
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
    lattice_list.push_back(lattice);
  }

  RotateLatticeVector(lattice_list,
                      GetLatticePairRotationMatrix(config, lattice_id_jump_pair));
  std::sort(lattice_list.begin(), lattice_list.end(),
            [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
              return LatticeSortCompare(lhs, rhs);
            });

  return lattice_list;
}

static std::unordered_map<
    cfg::ElementCluster, int, boost::hash<cfg::ElementCluster> > InitializeClusterHashMap(
    const std::set<Element> &type_set) {
  std::unordered_map<cfg::ElementCluster, int,
                     boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap;

  for (const auto &element1: type_set) {
    initialized_cluster_hashmap[cfg::ElementCluster(0, element1)] = 0;
    for (const auto &element2: type_set) {
      if (element2 == ElementName::X) {
        continue;
      }
      if (element1 == ElementName::X && element2.GetString()[0] == 'p') {
        continue;
      }
      for (size_t label = 1; label <= 3; ++label) {
        initialized_cluster_hashmap[cfg::ElementCluster(label, element1, element2)] = 0;
      }
      for (const auto &element3: type_set) {
        if (element3 == ElementName::X || element3.GetString()[0] == 'p') {
          continue;
        }
        for (size_t label = 4; label < 8; ++label) {
          initialized_cluster_hashmap[cfg::ElementCluster(label, element1, element2, element3)] = 0;
        }
      }
    }
  }
  return initialized_cluster_hashmap;
}

static int GetLabel(const std::vector<size_t> &lattice_index_list, const cfg::Config &config) {
  if (lattice_index_list.size() == 1) {
    return 0;
  }
  if (lattice_index_list.size() == 2) {
    return cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[0],
                                                 lattice_index_list[1],
                                                 config);
  }
  if (lattice_index_list.size() == 3) {
    std::vector<int>
        bond_label_list{cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[0],
                                                              lattice_index_list[1],
                                                              config),
                        cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[1],
                                                              lattice_index_list[2],
                                                              config),
                        cfg::FindDistanceLabelBetween_Lattice(lattice_index_list[2],
                                                              lattice_index_list[0],
                                                              config)};
    std::sort(bond_label_list.begin(), bond_label_list.end());
    if (bond_label_list[0] == 1 && bond_label_list[1] == 1 && bond_label_list[2] == 1) {
      return 4;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 1 && bond_label_list[2] == 2) {
      return 5;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 1 && bond_label_list[2] == 3) {
      return 6;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 2 && bond_label_list[2] == 3) {
      return 7;
    } else if (bond_label_list[0] == 1 && bond_label_list[1] == 3 && bond_label_list[2] == 3) {
      return 8;
    } else if (bond_label_list[0] == 2 && bond_label_list[1] == 3 && bond_label_list[2] == 3) {
      return 9;
    } else if (bond_label_list[0] == 3 && bond_label_list[1] == 3 && bond_label_list[2] == 3) {
      return 10;
    }
  }
  return -1;
}

template<size_t DataSize>
static void GetParametersMappingFromLatticeClusterVectorHelper(
    std::vector<cfg::LatticeCluster<DataSize>> &&cluster_vector,
    std::vector<std::vector<std::vector<size_t> > > &cluster_mapping
) {
  std::vector<std::vector<size_t> > cluster_index_vector;
  for (const auto &cluster: cluster_vector) {
    auto cluster_index = cluster.GetIndexVector();
    cluster_index_vector.push_back(cluster_index);
  }
  cluster_mapping.push_back(cluster_index_vector);
}

static std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingState(
    const cfg::Config &config) {
  const std::pair<size_t, size_t>
      lattice_id_jump_pair = {0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto lattice_vector = GetSortedLatticeVectorState(config, lattice_id_jump_pair);
  std::vector<std::vector<std::vector<size_t> > > cluster_mapping{};
  std::vector<Singlet_State_t> singlet_vector;
  std::vector<Pair_State_t> first_pair_vector;
  std::vector<Pair_State_t> second_pair_vector;
  std::vector<Pair_State_t> third_pair_vector;
  std::vector<Triplet_State_t> first_first_first_triplets_vector;
  std::vector<Triplet_State_t> first_first_second_triplets_vector;
  std::vector<Triplet_State_t> first_first_third_triplets_vector;
  std::vector<Triplet_State_t> first_second_third_triplets_vector;
  std::vector<Triplet_State_t> first_third_third_triplets_vector;
  std::vector<Triplet_State_t> second_third_third_triplets_vector;
  std::vector<Triplet_State_t> third_third_third_triplets_vector;

  for (size_t index1 = 0; index1 < lattice_vector.size(); ++index1) {
    cfg::Lattice lattice1(lattice_vector[index1]);
    const size_t id1 = lattice1.GetId();
    lattice1.SetId(index1);
    if (id1 == lattice_id_jump_pair.first || id1 == lattice_id_jump_pair.second) {
      singlet_vector.emplace_back(std::array<cfg::Lattice, 1>{lattice1});
    }
    for (size_t index2 = 0; index2 < index1; ++index2) {
      cfg::Lattice lattice2(lattice_vector[index2]);
      const size_t id2 = lattice2.GetId();
      lattice2.SetId(index2);
      if (id1 == lattice_id_jump_pair.first || id1 == lattice_id_jump_pair.second
          || id2 == lattice_id_jump_pair.first || id2 == lattice_id_jump_pair.second) {
        switch (GetLabel({id1, id2}, config)) {
          case 1: first_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 2: second_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          case 3: third_pair_vector.emplace_back(std::array<cfg::Lattice, 2>{lattice1, lattice2});
            break;
          default:continue;
        }
      }
      for (size_t index3 = 0; index3 < index2; ++index3) {
        cfg::Lattice lattice3(lattice_vector[index3]);
        const size_t id3 = lattice3.GetId();
        lattice3.SetId(index3);
        if (id1 == lattice_id_jump_pair.first || id1 == lattice_id_jump_pair.second
            || id2 == lattice_id_jump_pair.first || id2 == lattice_id_jump_pair.second
            || id3 == lattice_id_jump_pair.first || id3 == lattice_id_jump_pair.second) {
          switch (GetLabel({id1, id2, id3}, config)) {
            case 4:
              first_first_first_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 5:
              first_first_second_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 6:
              first_first_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 7:
              first_second_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 8:
              first_third_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 9:
              second_third_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
            case 10:
              third_third_third_triplets_vector.emplace_back(
                  std::array<cfg::Lattice, 3>{lattice1, lattice2, lattice3});
              break;
          }
        }
      }
    }
  }
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(singlet_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(second_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(third_pair_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_first_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_second_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_first_third_triplets_vector), cluster_mapping);
  GetParametersMappingFromLatticeClusterVectorHelper(
      std::move(first_second_third_triplets_vector), cluster_mapping);
  // GetParametersMappingFromLatticeClusterVectorHelper(
  //     std::move(first_third_third_triplets_vector), cluster_mapping);
  // GetParametersMappingFromLatticeClusterVectorHelper(
  //     std::move(second_third_third_triplets_vector), cluster_mapping);
  // GetParametersMappingFromLatticeClusterVectorHelper(
  //     std::move(third_third_third_triplets_vector), cluster_mapping);
  return cluster_mapping;
}

EnergyPredictor::EnergyPredictor(const std::string &predictor_filename,
                                 const cfg::Config &reference_config,
                                 const std::set<Element> &type_set)
    : cluster_mapping_(GetClusterParametersMappingState(reference_config)) {
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
    type_set_copy.emplace(ElementName::X);
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
      site_bond_cluster_hashmap_[{i, j}] = lattice_id_vector;
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
  auto migration_element = config.GetElementAtLatticeId(lattice_id_jump_pair.second);
  const auto &lattice_id_vector = site_bond_cluster_hashmap_.at(lattice_id_jump_pair);
  const auto &initialized_cluster_hashmap
      = element_initialized_cluster_hashmap_.at(migration_element);

  auto start_hashmap(initialized_cluster_hashmap);
  auto end_hashmap(initialized_cluster_hashmap);
  auto transition_hashmap(initialized_cluster_hashmap);
  size_t label = 0;
  for (const auto &cluster_vector: cluster_mapping_) {
    for (const auto &cluster: cluster_vector) {
      std::vector<Element> element_vector_start, element_vector_end, element_vector_transition;
      element_vector_start.reserve(cluster.size());
      element_vector_end.reserve(cluster.size());
      element_vector_transition.reserve(cluster.size());
      for (auto index: cluster) {
        size_t lattice_id = lattice_id_vector[index];
        element_vector_start.push_back(config.GetElementAtLatticeId(lattice_id));
        if (lattice_id == lattice_id_jump_pair.first) {
          element_vector_end.push_back(migration_element);
          element_vector_transition.push_back(migration_element.GetPseudo());
          continue;
        } else if (lattice_id == lattice_id_jump_pair.second) {
          element_vector_end.emplace_back(ElementName::X);
          element_vector_transition.push_back(migration_element.GetPseudo());
          continue;
        }
        element_vector_end.push_back(config.GetElementAtLatticeId(lattice_id));
        element_vector_transition.push_back(config.GetElementAtLatticeId(lattice_id));
      }
      start_hashmap[cfg::ElementCluster(label, element_vector_start)]++;
      end_hashmap[cfg::ElementCluster(label, element_vector_end)]++;
      transition_hashmap[cfg::ElementCluster(label, element_vector_transition)]++;
    }
    label++;
  }

  std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap.begin(), initialized_cluster_hashmap.end());
  std::vector<double> de_encode, e0_encode;
  de_encode.reserve(ordered.size());
  e0_encode.reserve(ordered.size());
  for (const auto &cluster_count: ordered) {
    const auto &cluster = cluster_count.first;
    auto start = static_cast<double>(start_hashmap.at(cluster));
    auto end = static_cast<double>(end_hashmap.at(cluster));
    auto transition = static_cast<double>(transition_hashmap.at(cluster));
    double total_bond{};
    switch (cluster.GetLabel()) {
      case 0:total_bond = 256;
        break;
      case 1: total_bond = 1536;
        break;
      case 2: total_bond = 768;
        break;
      case 3: total_bond = 3072;
        break;
      case 4: total_bond = 2048;
        break;
      case 5: total_bond = 3072;
        break;
      case 6: total_bond = 6144;
        break;
      case 7: total_bond = 6144;
        break;
      case 8: total_bond = 6144;
        break;
      case 9: total_bond = 6144;
        break;
      case 10: total_bond = 2048;
        break;
    }
    de_encode.push_back((end - start) / total_bond);
    e0_encode.push_back((transition - 0.5 * (end + start)) / total_bond);
  }

  const auto &theta_element = element_theta_.at(migration_element);
  double e0 = 0, dE = 0;

  const size_t cluster_size = theta_element.size();
  for (size_t i = 0; i < cluster_size; ++i) {
    dE += base_theta_[i] * de_encode[i];
    e0 += theta_element[i] * e0_encode[i];
  }

  return {e0 + dE / 2, dE};
}
} // namespace pred
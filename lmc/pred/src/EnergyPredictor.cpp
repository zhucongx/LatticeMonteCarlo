#include "EnergyPredictor.h"
#include <omp.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <ranges>

namespace pred {
EnergyPredictor::EnergyPredictor(const std::string &predictor_filename,
                                 std::set<Element> element_set)
    : element_set_(std::move(element_set)) {
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  initialized_cluster_hashmap_ = InitializeClusterHashMap(element_set_copy);

  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());

  std::ifstream ifs(predictor_filename, std::ifstream::in);
  if (!ifs) {
    throw std::runtime_error("Cannot open " + predictor_filename);
  }
  nlohmann::json all_parameters;
  ifs >> all_parameters;

  for (const auto &[element, parameters] : all_parameters.items()) {
    if (element == "Base") {
      base_theta_ = JsonToEigenVector(parameters.at("theta"));
    }
  }
}
EnergyPredictor::~EnergyPredictor() = default;
std::vector<double> EnergyPredictor::GetEncode(const cfg::Config &config) const {
  auto cluster_hashmap(initialized_cluster_hashmap_);
  static const std::vector<size_t>
      cluster_counter{256, 3072, 1536, 6144, 12288, 6144, 12288, 6144, 12288, 12288, 12288};
  for (size_t atom_id1 = 0; atom_id1 < config.GetNumAtoms(); ++atom_id1) {
    const size_t lattice_id1 = config.GetLatticeIdFromAtomId(atom_id1);
    Element element1 = config.GetAtomVector()[atom_id1].GetElement();
    cluster_hashmap[cfg::ElementCluster(0, element1)]++;
    for (size_t lattice_id2 : config.GetFirstNeighborsAdjacencyList()[lattice_id1]) {
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(1, element1, element2)]++;
      for (size_t lattice_id3 : config.GetFirstNeighborsAdjacencyList()[lattice_id2]) {
        Element element3 = config.GetElementAtLatticeId(lattice_id3);
        if (std::ranges::find(config.GetFirstNeighborsAdjacencyList()[lattice_id1],
                      lattice_id3)
            != config.GetFirstNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(4, element1, element2, element3)]++;
        } else if (std::ranges::find(config.GetSecondNeighborsAdjacencyList()[lattice_id1],
                             lattice_id3)
            != config.GetSecondNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(5, element1, element2, element3)]++;
        } else if (std::ranges::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1],
                             lattice_id3)
            != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(6, element1, element2, element3)]++;
        }
      }
      for (size_t lattice_id3 : config.GetSecondNeighborsAdjacencyList()[lattice_id2]) {
        Element element3 = config.GetElementAtLatticeId(lattice_id3);
        if (std::ranges::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1],
                      lattice_id3)
            != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(7, element1, element2, element3)]++;
        }
      }
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(8, element1, element2, element3)]++;
      //   }
      // }
    }
    for (size_t lattice_id2 : config.GetSecondNeighborsAdjacencyList()[lattice_id1]) {
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(2, element1, element2)]++;
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(9, element1, element2, element3)]++;
      //   }
      // }
    }
    for (size_t lattice_id2 : config.GetThirdNeighborsAdjacencyList()[lattice_id1]) {
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(3, element1, element2)]++;
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(10, element1, element2, element3)]++;
      //   }
      // }
    }
  }

  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  std::vector<double> energy_encode;
  for (const auto &cluster_count : ordered) {
    const auto &cluster = cluster_count.first;
    auto count_bond = static_cast<double>(cluster_hashmap.at(cluster));
    auto total_bond = static_cast<double>(cluster_counter[static_cast<size_t>(cluster.GetLabel())]);
    energy_encode.push_back(count_bond / total_bond);
  }
  return energy_encode;
}
std::vector<double> EnergyPredictor::GetEncodeOfCluster(
    const cfg::Config &config, const std::vector<size_t> &atom_id_list) const {
  auto cluster_hashmap(initialized_cluster_hashmap_);
  std::unordered_set<size_t> lattice_id_hashset;
  for (auto atom_id : atom_id_list) {
    auto lattice_id = config.GetLatticeIdFromAtomId(atom_id);
    lattice_id_hashset.insert(lattice_id);
    for (auto neighbor_lattice_id : config.GetFirstNeighborsAdjacencyList().at(lattice_id)) {
      lattice_id_hashset.insert(neighbor_lattice_id);
    }
    for (auto neighbor_lattice_id : config.GetSecondNeighborsAdjacencyList().at(lattice_id)) {
      lattice_id_hashset.insert(neighbor_lattice_id);
    }
    for (auto neighbor_lattice_id : config.GetThirdNeighborsAdjacencyList().at(lattice_id)) {
      lattice_id_hashset.insert(neighbor_lattice_id);
    }
  }
  const std::vector<size_t>
      cluster_counter{256, 3072, 1536, 6144, 12288, 6144, 12288, 6144, 12288, 12288, 12288};
  for (size_t lattice_id1 : lattice_id_hashset) {
    Element element1 = config.GetElementAtLatticeId(lattice_id1);
    cluster_hashmap[cfg::ElementCluster(0, element1)]++;
    for (size_t lattice_id2 : config.GetFirstNeighborsAdjacencyList()[lattice_id1]) {
      if (lattice_id_hashset.find(lattice_id2) == lattice_id_hashset.end()) { continue; }
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(1, element1, element2)]++;
      for (size_t lattice_id3 : config.GetFirstNeighborsAdjacencyList()[lattice_id2]) {
        if (!lattice_id_hashset.contains(lattice_id3)) { continue; }
        Element element3 = config.GetElementAtLatticeId(lattice_id3);
        if (std::ranges::find(config.GetFirstNeighborsAdjacencyList()[lattice_id1],
                      lattice_id3)
            != config.GetFirstNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(4, element1, element2, element3)]++;
        } else if (std::ranges::find(config.GetSecondNeighborsAdjacencyList()[lattice_id1],
                             lattice_id3)
            != config.GetSecondNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(5, element1, element2, element3)]++;
        } else if (std::ranges::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1],
                             lattice_id3)
            != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(6, element1, element2, element3)]++;
        }
      }
      for (size_t lattice_id3 : config.GetSecondNeighborsAdjacencyList()[lattice_id2]) {
        if (!lattice_id_hashset.contains(lattice_id3)) { continue; }
        Element element3 = config.GetElementAtLatticeId(lattice_id3);
        if (std::ranges::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1],
                      lattice_id3)
            != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
          cluster_hashmap[cfg::ElementCluster(7, element1, element2, element3)]++;
        }
      }
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   if (lattice_id_set.find(lattice_id3) == lattice_id_set.end()) { continue; }
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(8, element1, element2, element3)]++;
      //   }
      // }
    }
    for (size_t lattice_id2 : config.GetSecondNeighborsAdjacencyList()[lattice_id1]) {
      if (!lattice_id_hashset.contains(lattice_id2)) { continue; }
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(2, element1, element2)]++;
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   if (lattice_id_set.find(lattice_id3) == lattice_id_set.end()) { continue; }
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(9, element1, element2, element3)]++;
      //   }
      // }
    }
    for (size_t lattice_id2 : config.GetThirdNeighborsAdjacencyList()[lattice_id1]) {
      if (!lattice_id_hashset.contains(lattice_id2)) { continue; }
      Element element2 = config.GetElementAtLatticeId(lattice_id2);
      cluster_hashmap[cfg::ElementCluster(3, element1, element2)]++;
      // for (size_t lattice_id3: config.GetThirdNeighborsAdjacencyList()[lattice_id2]) {
      //   if (lattice_id_set.find(lattice_id3) == lattice_id_set.end()) { continue; }
      //   Element element3 = config.GetElementAtLatticeId(lattice_id3);
      //   if (std::find(config.GetThirdNeighborsAdjacencyList()[lattice_id1].begin(),
      //                 config.GetThirdNeighborsAdjacencyList()[lattice_id1].end(),
      //                 lattice_id3)
      //       != config.GetThirdNeighborsAdjacencyList()[lattice_id1].end()) {
      //     cluster_hashmap[cfg::ElementCluster(10, element1, element2, element3)]++;
      //   }
      // }
    }
  }

  const std::map<cfg::ElementCluster, int>
      ordered(initialized_cluster_hashmap_.begin(), initialized_cluster_hashmap_.end());
  std::vector<double> energy_encode;
  for (const auto &cluster: ordered | std::views::keys) {
    auto count_bond = static_cast<double>(cluster_hashmap.at(cluster));
    auto total_bond = static_cast<double>(cluster_counter[static_cast<size_t>(cluster.GetLabel())]);
    energy_encode.push_back(count_bond / total_bond);
  }
  return energy_encode;
}
double EnergyPredictor::GetEnergy(const cfg::Config &config) const {
  const auto encode = GetEncode(config);
  const Eigen::Map<const Eigen::VectorXd> encode_vec(encode.data(), static_cast<Eigen::Index>(encode.size()));
  return base_theta_.dot(encode_vec);
}
double EnergyPredictor::GetEnergyOfCluster(
    const cfg::Config &config,
    const std::vector<size_t> &atom_id_list) const {
  const auto encode = GetEncodeOfCluster(config, atom_id_list);
  const Eigen::Map<const Eigen::VectorXd> encode_vec(encode.data(), static_cast<Eigen::Index>(encode.size()));
  return base_theta_.dot(encode_vec);
}
// Todo: add a function to get the chemical potential of a given composition
// std::map<Element, double> EnergyPredictor::GetChemicalPotential() const {
//   std::map<Element, double> chemical_potential;
//   for (auto element: element_set_) {
//     auto config = cfg::GenerateFCC({10, 10, 10}, element);
//     chemical_potential[element] = GetEnergy(config) / static_cast<double>(config.GetNumAtoms());
//     std::cerr << element.GetString() << " " << chemical_potential[element] << std::endl;
//   }
//   return chemical_potential;
// }

std::map<Element, double> EnergyPredictor::GetChemicalPotential(Element solvent_element) const {
  std::map<Element, double> chemical_potential;
  const double solvent_energy = GetEnergy(cfg::GenerateFCC({15, 15, 15}, solvent_element));
  chemical_potential[solvent_element] = 0;
  auto element_set_copy(element_set_);
  element_set_copy.emplace(ElementName::X);
  for (auto element : element_set_copy) {
    if (element == solvent_element) {
      continue;
    }
    auto config = cfg::GenerateFCC({15, 15, 15}, solvent_element);
    config.SetAtomElementTypeAtAtom(0, element);
    chemical_potential[element] = GetEnergy(config) - solvent_energy;
  }
  // for (auto [element, potential] : chemical_potential) {
  //   std::cout << "Chemical potential: " << element.GetString() << " " << potential << std::endl;
  // }
  return chemical_potential;
}
} // pred

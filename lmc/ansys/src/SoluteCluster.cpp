#include "SoluteCluster.h"

#include "ShortRangeOrder.h"

#include <Eigen/Dense>
#include <cmath>
#include <filesystem>
#include <omp.h>
#include <queue>
#include <unordered_map>
#include <utility>

namespace ansys {

// template<typename T>
// static void AppendInfoToAuxiliaryListsRepeat(std::map<std::string, cfg::Config::VectorVariant> &auxiliary_lists,
//                                              const std::string &key,
//                                              T value,
//                                              size_t repeat)
// {
//   const auto it = auxiliary_lists.find(key);
//   if (it == auxiliary_lists.end()) {
//     auxiliary_lists[key] = std::vector<T>();
//   }
//   std::visit(
//       [&value, &repeat](auto &&vec) {
//         using VecType = std::decay_t<decltype(vec)>;
//         if constexpr (std::is_same_v<VecType, std::vector<T>>) {
//           for (size_t i = 0; i < repeat; ++i) {
//             vec.push_back(value);
//           }
//         } else {
//           throw std::runtime_error("Type mismatch in auxiliary_lists");
//         }
//       },
//       it->second);
// }

// template<typename T>
// void ModifyElement(std::map<std::string, cfg::Config::VectorVariant> &auxiliary_lists,
//                    const std::string &key,
//                    size_t index,
//                    const T &new_value)
// {
//   const auto it = auxiliary_lists.find(key);
//   if (it == auxiliary_lists.end()) {
//     throw std::runtime_error("Key not found in auxiliary_lists");
//   }
//   std::visit(
//       [index, &new_value](auto &vec) {
//         using VecType = std::decay_t<decltype(vec)>;
//         if constexpr (std::is_same_v<VecType, std::vector<T>>) {
//           if (index < vec.size()) {
//             vec[index] = new_value;
//           }
//         } else {
//           throw std::runtime_error("Type mismatch in auxiliary_lists");
//         }
//       },
//       it->second);
// }


SoluteCluster::SoluteCluster(const cfg::Config &config,
                             Element solvent_atom_type,
                             std::set<Element> element_set,
                             size_t smallest_cluster_criteria,
                             size_t solvent_bond_criteria,
                             const pred::EnergyPredictor &energy_predictor,
                             const std::map<Element, double> &chemical_potential_map)
    : config_(config),
      solvent_config_(config),
      solvent_element_(solvent_atom_type),
      element_set_(std::move(element_set)),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      solvent_bond_criteria_(solvent_bond_criteria),
      energy_predictor_(energy_predictor),
      chemical_potential_map_(chemical_potential_map) {
  for (size_t atom_id = 0; atom_id < solvent_config_.GetNumAtoms(); ++atom_id) {
    solvent_config_.SetAtomElementTypeAtAtom(atom_id, solvent_atom_type);
  }
}

template<typename T>
void ModifyElementFromVector(std::map<std::string, cfg::Config::VectorVariant> &auxiliary_lists,
                             const std::string &key,
                             const std::vector<size_t> &indices,
                             const T &new_value,
                             const std::vector<T> default_vector) {
  if (auxiliary_lists.find(key) == auxiliary_lists.end()) {
    auxiliary_lists[key] = default_vector;
  }
  std::visit(
      [&indices, &new_value](auto &vec) {
        using VecType = std::decay_t<decltype(vec)>;
        if constexpr (std::is_same_v<VecType, std::vector<T>>) {
          for (const auto index: indices) {
            if (index < vec.size()) {
              vec[index] = new_value;
            }
          }
        } else {
          throw std::runtime_error("Type mismatch in auxiliary_lists");
        }
      },
      auxiliary_lists[key]);
}

std::pair<nlohmann::json, std::map<std::string, cfg::Config::VectorVariant>> SoluteCluster::GetClustersInfo() {
  const auto cluster_to_atom_vector = FindAtomListOfClusters();

  nlohmann::json clusters_info_array = nlohmann::json::array();

  std::map<std::string, cfg::Config::VectorVariant> auxiliary_lists{};
  for (const auto &element: element_set_) {
    auxiliary_lists["cluster_" + element.GetString()] = std::vector<size_t>(config_.GetNumAtoms(), 0);
  }
  for (size_t id = 0; id < cluster_to_atom_vector.size(); ++id) {
    const auto &cluster_atom_id_list = cluster_to_atom_vector[id];

    nlohmann::json cluster_info = nlohmann::json::object();

    const size_t size = cluster_atom_id_list.size();
    const auto element_number = GetElementsNumber(cluster_atom_id_list);
    const size_t size_without_vacancy = size - element_number.at("X");

    cluster_info["cluster_atom_id_list"] = cluster_atom_id_list;

    cluster_info["cluster_id"] = id + 1;

    ModifyElementFromVector(
        auxiliary_lists, "cluster_id", cluster_atom_id_list, id + 1, std::vector<size_t>(config_.GetNumAtoms(), 0));

    cluster_info["cluster_size"] = size_without_vacancy;
    ModifyElementFromVector(auxiliary_lists,
                            "cluster_size",
                            cluster_atom_id_list,
                            size_without_vacancy,
                            std::vector<size_t>(config_.GetNumAtoms(), 0));

    const double total_volume = static_cast<double>(size) * std::pow(constants::kLatticeConstant, 3) /
        static_cast<double>(constants::kNumAtomsPerCell);
    const double effective_radius = std::pow(3 * total_volume / (4 * M_PI), 1.0 / 3.0);
    cluster_info["effective_radius"] = effective_radius;
    ModifyElementFromVector(auxiliary_lists,
                            "effective_radius",
                            cluster_atom_id_list,
                            effective_radius,
                            std::vector<double>(config_.GetNumAtoms(), NAN));

    cluster_info["elements_number"] = element_number;
    for (const auto &[element, number]: element_number) {
      ModifyElementFromVector(auxiliary_lists,
                              "cluster_" + element,
                              cluster_atom_id_list,
                              number,
                              std::vector<size_t>(config_.GetNumAtoms(), 0));
    }

    const double cluster_mass = GetMass(cluster_atom_id_list);
    cluster_info["cluster_mass"] = cluster_mass;
    ModifyElementFromVector(auxiliary_lists,
                            "cluster_mass",
                            cluster_atom_id_list,
                            cluster_mass,
                            std::vector<double>(config_.GetNumAtoms(), NAN));

    const auto cluster_energy = GetFormationEnergy(cluster_atom_id_list) / static_cast<double>(size_without_vacancy);
    cluster_info["cluster_energy"] = cluster_energy;
    ModifyElementFromVector(auxiliary_lists,
                            "cluster_energy",
                            cluster_atom_id_list,
                            cluster_energy,
                            std::vector<double>(config_.GetNumAtoms(), NAN));

    const auto geometry_center = GetGeometryCenter(cluster_atom_id_list);
    cluster_info["geometry_center"] = geometry_center;
    ModifyElementFromVector(auxiliary_lists,
                            "geometry_center",
                            cluster_atom_id_list,
                            geometry_center,
                            std::vector<Vector_t>(config_.GetNumAtoms(), {NAN, NAN, NAN}));

    const auto mass_center = GetMassCenter(cluster_atom_id_list);
    cluster_info["mass_center"] = mass_center;
    ModifyElementFromVector(auxiliary_lists,
                            "mass_center",
                            cluster_atom_id_list,
                            mass_center,
                            std::vector<Vector_t>(config_.GetNumAtoms(), {NAN, NAN, NAN}));

    const auto mass_gyration_tensor = GetMassGyrationTensor(cluster_atom_id_list, mass_center);
    cluster_info["mass_gyration_tensor"] = mass_gyration_tensor;
    // ModifyElementFromVector(auxiliary_lists, "mass_gyration_tensor", cluster_atom_id_list, mass_gyration_tensor);

    Eigen::Matrix3d mass_gyration_tensor_eigen{
        {mass_gyration_tensor[0][0], mass_gyration_tensor[0][1], mass_gyration_tensor[0][2]},
        {mass_gyration_tensor[1][0], mass_gyration_tensor[1][1], mass_gyration_tensor[1][2]},
        {mass_gyration_tensor[2][0], mass_gyration_tensor[2][1], mass_gyration_tensor[2][2]}};
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(mass_gyration_tensor_eigen);
    if (eigen_solver.info() != Eigen::Success) {
      throw std::runtime_error("Eigen solver failed");
    }
    const auto &eigenvalues = eigen_solver.eigenvalues();
    const auto mass_gyration_radius = std::sqrt(eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
    cluster_info["mass_gyration_radius"] = mass_gyration_radius;
    ModifyElementFromVector(auxiliary_lists,
                            "mass_gyration_radius",
                            cluster_atom_id_list,
                            mass_gyration_radius,
                            std::vector<double>(config_.GetNumAtoms(), NAN));

    const auto asphericity = eigenvalues[2] - 0.5 * (eigenvalues[0] + eigenvalues[1]);
    cluster_info["asphericity"] = asphericity;
    ModifyElementFromVector(auxiliary_lists,
                            "cluster_asphericity",
                            cluster_atom_id_list,
                            asphericity,
                            std::vector<double>(config_.GetNumAtoms(), NAN));

    const auto acylindricity = eigenvalues[1] - eigenvalues[0];
    cluster_info["acylindricity"] = acylindricity;
    ModifyElementFromVector(auxiliary_lists,
                            "cluster_acylindricity",
                            cluster_atom_id_list,
                            acylindricity,
                            std::vector<double>(config_.GetNumAtoms(), NAN));

    const auto anisotropy = 1.5 *
            (std::pow(eigenvalues[0], 2) + std::pow(eigenvalues[1], 2) + std::pow(eigenvalues[2], 2)) /
            std::pow(eigenvalues[0] + eigenvalues[1] + eigenvalues[2], 2) -
        0.5;
    cluster_info["anisotropy"] = anisotropy;
    ModifyElementFromVector(auxiliary_lists,
                            "cluster_anisotropy",
                            cluster_atom_id_list,
                            anisotropy,
                            std::vector<double>(config_.GetNumAtoms(), NAN));

    cluster_info["mass_inertia_tensor"] = GetMassInertiaTensor(cluster_atom_id_list, mass_center);

    clusters_info_array.push_back(cluster_info);
  }

  return std::make_pair(clusters_info_array, auxiliary_lists);
}

std::unordered_set<size_t> SoluteCluster::FindSoluteAtomIndexes() const {
  std::unordered_set<size_t> solute_atoms_hashset;
  for (const auto &atom: config_.GetAtomVector()) {
    if (atom.GetElement() == solvent_element_) {
      continue;
    }
    solute_atoms_hashset.insert(atom.GetId());
  }
  return solute_atoms_hashset;
}

std::vector<std::vector<size_t>>
SoluteCluster::FindAtomListOfClustersBFSHelper(std::unordered_set<size_t> unvisited_atoms_id_set) const {
  std::vector<std::vector<size_t>> cluster_atom_list;
  std::queue<size_t> visit_id_queue;
  size_t atom_id;

  std::unordered_set<size_t>::iterator it;
  while (!unvisited_atoms_id_set.empty()) {
    // Find next element
    it = unvisited_atoms_id_set.begin();
    visit_id_queue.push(*it);
    unvisited_atoms_id_set.erase(it);

    std::vector<size_t> atom_list_of_one_cluster;
    while (!visit_id_queue.empty()) {
      atom_id = visit_id_queue.front();
      visit_id_queue.pop();

      atom_list_of_one_cluster.push_back(atom_id);
      for (const auto &neighbors_list: {
               config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id),
               // config_.GetSecondNeighborsAtomIdVectorOfAtom(atom_id)
           }) {
        for (auto neighbor_id: neighbors_list) {
          it = unvisited_atoms_id_set.find(neighbor_id);
          if (it != unvisited_atoms_id_set.end()) {
            visit_id_queue.push(*it);
            unvisited_atoms_id_set.erase(it);
          }
        }
      }
    }
    cluster_atom_list.push_back(atom_list_of_one_cluster);
  }
  return cluster_atom_list;
}

std::vector<std::vector<size_t>> SoluteCluster::FindAtomListOfClusters() const {
  auto cluster_atom_list = FindAtomListOfClustersBFSHelper(FindSoluteAtomIndexes());

  // remove small clusters
  auto it = cluster_atom_list.begin();
  while (it != cluster_atom_list.end()) {
    if (it->size() < smallest_cluster_criteria_) {
      it = cluster_atom_list.erase(it);
    } else {
      ++it;
    }
  }
  // add solvent neighbors
  for (auto &cluster: cluster_atom_list) {
    std::unordered_set<size_t> cluster_set(cluster.begin(), cluster.end());
    std::unordered_set<size_t> neighbor_set;
    for (auto atom_id: cluster) {
      neighbor_set.insert(atom_id);
      for (auto neighbor_id: config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id)) {
        neighbor_set.insert(neighbor_id);
      }
    }

    for (auto atom_id: neighbor_set) {
      size_t neighbor_bond_count = 0;
      for (auto neighbor_id: config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id)) {
        if (cluster_set.find(neighbor_id) != cluster_set.end() &&
            config_.GetAtomVector()[neighbor_id].GetElement() != solvent_element_) {
          neighbor_bond_count++;
        }
      }
      if (neighbor_bond_count > solvent_bond_criteria_) {
        cluster_set.insert(atom_id);
      }
    }
    cluster = std::vector<size_t>(cluster_set.begin(), cluster_set.end());
  }

  // sort by size
  std::sort(cluster_atom_list.begin(),
            cluster_atom_list.end(),
            [](const std::vector<size_t> &a, const std::vector<size_t> &b) {
              return a.size() > b.size();
            });

  return cluster_atom_list;
}

std::map<std::string, size_t> SoluteCluster::GetElementsNumber(const std::vector<size_t> &cluster_atom_id_list) const {
  // initialize map with all the element, because some cluster may not have all types of element
  std::map<std::string, size_t> num_atom_in_one_cluster{{"X", 0}};
  for (const auto &element: element_set_) {
    num_atom_in_one_cluster[element.GetString()] = 0;
  }
  for (const auto &atom_id: cluster_atom_id_list) {
    num_atom_in_one_cluster.at(config_.GetAtomVector()[atom_id].GetElement().GetString())++;
  }
  return num_atom_in_one_cluster;
}

double SoluteCluster::GetMass(const std::vector<size_t> &cluster_atom_id_list) const {
  double sum_mass = 0;
  for (const auto &atom_id: cluster_atom_id_list) {
    sum_mass += config_.GetAtomVector()[atom_id].GetElement().GetMass();
  }
  return sum_mass;
}

double SoluteCluster::GetFormationEnergy(const std::vector<size_t> &cluster_atom_id_list) const {
  cfg::Config solute_config(solvent_config_);
  double energy_change_solution_to_pure_solvent = 0;
  for (size_t atom_id: cluster_atom_id_list) {
    Element element = config_.GetElementAtAtomId(atom_id);
    solute_config.SetAtomElementTypeAtAtom(atom_id, element);
    energy_change_solution_to_pure_solvent += chemical_potential_map_.at(element);
  }
  const double energy_change_cluster_to_pure_solvent =
      energy_predictor_.GetEnergyOfCluster(solute_config, cluster_atom_id_list) -
      energy_predictor_.GetEnergyOfCluster(solvent_config_, cluster_atom_id_list);
  return energy_change_cluster_to_pure_solvent - energy_change_solution_to_pure_solvent;
}

Vector_t SoluteCluster::GetGeometryCenter(const std::vector<size_t> &cluster_atom_id_list) const {
  Vector_t geometry_center{};
  Vector_t sum_cos_theta{};
  Vector_t sum_sin_theta{};
  for (size_t atom_id: cluster_atom_id_list) {
    const auto relative_position =
        config_.GetLatticeVector()[config_.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    for (const auto kDim: All_Dimensions) {
      auto theta = relative_position[kDim] * 2 * M_PI;
      sum_cos_theta[kDim] += std::cos(theta);
      sum_sin_theta[kDim] += std::sin(theta);
    }
  }
  const auto cos_theta_bar = sum_cos_theta / static_cast<double>(cluster_atom_id_list.size());
  const auto sin_theta_bar = sum_sin_theta / static_cast<double>(cluster_atom_id_list.size());
  for (const auto kDim: All_Dimensions) {
    double theta_bar = std::atan2(-sin_theta_bar[kDim], -cos_theta_bar[kDim]) + M_PI;
    geometry_center[kDim] = theta_bar / (2 * M_PI);
  }
  return geometry_center * config_.GetBasis();
  ;    // Cartesian position
}

Vector_t SoluteCluster::GetMassCenter(const std::vector<size_t> &cluster_atom_id_list) const {
  // Cartesian position
  Vector_t mass_center{};
  Vector_t sum_cos_theta{};
  Vector_t sum_sin_theta{};
  double sum_mass = 0;
  for (size_t atom_id: cluster_atom_id_list) {
    auto relative_position = config_.GetLatticeVector()[config_.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    auto mass = config_.GetAtomVector()[atom_id].GetElement().GetMass();
    sum_mass += mass;
    for (const auto kDim: All_Dimensions) {
      auto theta = relative_position[kDim] * 2 * M_PI;
      sum_cos_theta[kDim] += std::cos(theta) * mass;
      sum_sin_theta[kDim] += std::sin(theta) * mass;
    }
  }
  auto cos_theta_bar = sum_cos_theta / sum_mass;
  auto sin_theta_bar = sum_sin_theta / sum_mass;
  for (const auto kDim: All_Dimensions) {
    double theta_bar = std::atan2(-sin_theta_bar[kDim], -cos_theta_bar[kDim]) + M_PI;
    mass_center[kDim] = theta_bar / (2 * M_PI);
  }
  return mass_center * config_.GetBasis();
}

Matrix_t SoluteCluster::GetMassGyrationTensor(const std::vector<size_t> &cluster_atom_id_list,
                                              const Vector_t &mass_center) const {
  const auto relative_mass_center = mass_center * InverseMatrix(config_.GetBasis());
  Matrix_t gyration_tensor{};
  double sum_mass = 0;
  for (size_t atom_id: cluster_atom_id_list) {
    const auto relative_position =
        config_.GetLatticeVector()[config_.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    auto mass = config_.GetAtomVector()[atom_id].GetElement().GetMass();
    sum_mass += mass;
    for (const auto kDim1: All_Dimensions) {
      auto r1 = relative_position[kDim1] - relative_mass_center[kDim1];
      while (r1 >= 0.5) {
        r1 -= 1;
      }
      while (r1 < -0.5) {
        r1 += 1;
      }
      for (const auto kDim2: All_Dimensions) {
        auto r2 = relative_position[kDim2] - relative_mass_center[kDim2];
        while (r2 >= 0.5) {
          r2 -= 1;
        }
        while (r2 < -0.5) {
          r2 += 1;
        }
        gyration_tensor[kDim1][kDim2] += r1 * r2 * mass;
      }
    }
  }
  gyration_tensor /= sum_mass;
  return {gyration_tensor[0] * config_.GetBasis() * config_.GetBasis(),
          gyration_tensor[1] * config_.GetBasis() * config_.GetBasis(),
          gyration_tensor[2] * config_.GetBasis() * config_.GetBasis()};
}

Matrix_t SoluteCluster::GetMassInertiaTensor(const std::vector<size_t> &cluster_atom_id_list,
                                             const Vector_t &mass_center) const {
  const auto relative_mass_center = mass_center * InverseMatrix(config_.GetBasis());
  Matrix_t inertia_tensor{};
  for (const size_t atom_id: cluster_atom_id_list) {
    const auto mass = config_.GetAtomVector()[atom_id].GetElement().GetMass();
    auto relative_distance = config_.GetLatticeVector()[config_.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition() -
        relative_mass_center;
    for (const auto kDim: All_Dimensions) {
      while (relative_distance[kDim] >= 0.5) {
        relative_distance[kDim] -= 1;
      }
      while (relative_distance[kDim] < -0.5) {
        relative_distance[kDim] += 1;
      }
    }
    auto cartesian_distance = relative_distance * config_.GetBasis();
    inertia_tensor[0][0] +=
        mass * (cartesian_distance[1] * cartesian_distance[1] + cartesian_distance[2] * cartesian_distance[2]);
    inertia_tensor[0][1] += -mass * (cartesian_distance[0] * cartesian_distance[1]);
    inertia_tensor[0][2] += -mass * (cartesian_distance[0] * cartesian_distance[2]);
    inertia_tensor[1][0] += -mass * (cartesian_distance[1] * cartesian_distance[0]);
    inertia_tensor[1][1] +=
        mass * (cartesian_distance[2] * cartesian_distance[2] + cartesian_distance[0] * cartesian_distance[0]);
    inertia_tensor[1][2] += -mass * (cartesian_distance[1] * cartesian_distance[2]);
    inertia_tensor[2][0] += -mass * (cartesian_distance[2] * cartesian_distance[0]);
    inertia_tensor[2][1] += -mass * (cartesian_distance[2] * cartesian_distance[1]);
    inertia_tensor[2][2] +=
        mass * (cartesian_distance[0] * cartesian_distance[0] + cartesian_distance[1] * cartesian_distance[1]);
  }
  return inertia_tensor;
}
}    // namespace ansys

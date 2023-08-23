/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 6/14/20 1:27 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/22/23 11:08 PM                                                          *
 **************************************************************************************************/

#include "SoluteCluster.h"

#include <queue>
#include <unordered_map>
#include <utility>
#include <filesystem>
#include <cmath>
#include <omp.h>
#include <boost/filesystem.hpp>
#include <Eigen/Dense>

using json = nlohmann::json;

template<typename T>
static void AppendInfoToAuxiliaryListsRepeat(std::map<std::string, Config::VectorVariant> &auxiliary_lists,
                                             const std::string &key,
                                             T value,
                                             size_t repeat) {
  if (auxiliary_lists.find(key) == auxiliary_lists.end()) {
    auxiliary_lists[key] = std::vector<T>();
  }
  std::visit([&value, &repeat](auto &&vec) {
    using VecType = std::decay_t<decltype(vec)>;
    if constexpr (std::is_same_v<VecType, std::vector<T>>) {
      for (size_t i = 0; i < repeat; ++i) {
        vec.push_back(value);
      }
    } else {
      throw std::runtime_error("Type mismatch in AppendInfoToAuxiliaryListsRepeat");
    }
  }, auxiliary_lists[key]);
}

SoluteCluster::SoluteCluster(const Config &config,
                             Element solvent_atom_type,
                             std::set<Element> element_set,
                             size_t smallest_cluster_criteria,
                             size_t solvent_bond_criteria,
                             const pred::EnergyPredictor &energy_estimator,
                             const std::map<Element, double> &chemical_potential_map)
    : config_(config),
      solvent_config_(config),
      solvent_element_(solvent_atom_type),
      element_set_(std::move(element_set)),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      solvent_bond_criteria_(solvent_bond_criteria),
      energy_estimator_(energy_estimator),
      chemical_potential_map_(chemical_potential_map) {
  for (size_t atom_id = 0; atom_id < solvent_config_.GetNumAtoms(); ++atom_id) {
    solvent_config_.SetElementOfAtom(atom_id, solvent_atom_type);
  }
}

json SoluteCluster::GetClustersInfoAndOutput(
    const std::string &output_folder,
    const std::string &output_name,
    const std::map<std::string, Config::ValueVariant> &global_info_map) {
  auto cluster_to_atom_vector = FindAtomListOfClusters();

  json clusters_info_array = json::array();
  std::vector<Eigen::Vector3d> relative_position_vector{};
  std::vector<Element> atom_vector{};
  for (auto &cluster_atom_id_list : cluster_to_atom_vector) {
    AppendAtomAndLatticeVector(cluster_atom_id_list, atom_vector, relative_position_vector);
    json cluster_info = json::object();

    const size_t size = cluster_atom_id_list.size();
    const auto element_number = GetElementsNumber(cluster_atom_id_list);
    const size_t size_without_vacancy = size - element_number.at("X");

    AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "cluster_size", size_without_vacancy, size);
    cluster_info["size"] = size_without_vacancy;

    cluster_info["elements_number"] = element_number;
    for (const auto &element_number_pair : element_number) {
      AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "cluster_" + element_number_pair.first,
                                       element_number_pair.second,
                                       size);
    }

    const double mass = GetMass(cluster_atom_id_list);
    cluster_info["mass"] = mass;

    const auto energy = GetFormationEnergy(cluster_atom_id_list) / static_cast<double> (size_without_vacancy);
    cluster_info["energy"] = energy;
    AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "cluster_energy", energy, size);

    const auto geometry_center = GetGeometryCenter(cluster_atom_id_list);
    cluster_info["geometry_center"] = geometry_center;

    const auto mass_center = GetMassCenter(cluster_atom_id_list);
    cluster_info["mass_center"] = mass_center;

    std::vector<std::vector<double> > matrix3d_helper(3, std::vector<double>(3, 0));
    Eigen::Map<Eigen::Matrix3d> mapped_data(matrix3d_helper[0].data());
    const auto mass_gyration_tensor = GetMassGyrationTensor(cluster_atom_id_list, mass_center);
    mapped_data = mass_gyration_tensor;
    cluster_info["mass_gyration_tensor"] = matrix3d_helper;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(mass_gyration_tensor);
    if (eigen_solver.info() != Eigen::Success) {
      throw std::runtime_error("Eigen solver failed");
    }
    const auto &eigenvalues = eigen_solver.eigenvalues();
    const auto mass_gyration_radius = std::sqrt(eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
    cluster_info["mass_gyration_radius"] = mass_gyration_radius;
    AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "cluster_mass_gyration_radius", mass_gyration_radius, size);

    const auto asphericity = eigenvalues[2] - 0.5 * (eigenvalues[0] + eigenvalues[1]);
    cluster_info["shape"]["asphericity"] = asphericity;
    AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "cluster_asphericity", asphericity, size);

    const auto acylindricity = eigenvalues[1] - eigenvalues[0];
    cluster_info["shape"]["acylindricity"] = acylindricity;
    AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "cluster_acylindricity", acylindricity, size);

    const auto anisotropy = 1.5
        * (std::pow(eigenvalues[0], 2) + std::pow(eigenvalues[1], 2) + std::pow(eigenvalues[2], 2))
        / std::pow(eigenvalues[0] + eigenvalues[1] + eigenvalues[2], 2) - 0.5;
    cluster_info["shape"]["anisotropy"] = anisotropy;
    AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "cluster_anisotropy", anisotropy, size);

    const auto mass_inertia_tensor = GetMassInertiaTensor(cluster_atom_id_list, mass_center);
    mapped_data = mass_inertia_tensor;
    cluster_info["mass_inertia_tensor"] = matrix3d_helper;

    // const auto first_pij = short_range_order.FindProbabilityCluster(1, cluster_atom_id_list);
    // cluster_info["pij"]["first"] = first_pij;
    // for (const auto &pair : first_pij) {
    //   AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "first_" + pair.first,
    //                                    static_cast<double> (pair.second),
    //                                    size);
    // }
    // const auto second_pij = short_range_order.FindProbabilityCluster(2, cluster_atom_id_list);
    // cluster_info["pij"]["second"] = second_pij;
    // for (const auto &pair : second_pij) {
    //   AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "second_" + pair.first,
    //                                    static_cast<double> (pair.second),
    //                                    size);
    // }
    // const auto third_pij = short_range_order.FindProbabilityCluster(3, cluster_atom_id_list);
    // cluster_info["pij"]["third"] = third_pij;
    // for (const auto &pair : third_pij) {
    //   AppendInfoToAuxiliaryListsRepeat(auxiliary_lists_, "third_" + pair.first,
    //                                    static_cast<double> (pair.second),
    //                                    size);
    // }

    clusters_info_array.push_back(cluster_info);
  }

  Eigen::Matrix3Xd relative_position_matrix =
      Eigen::Map<Eigen::Matrix3Xd>(reinterpret_cast<double *>(relative_position_vector.data()),
                                   3,
                                   static_cast<int>(relative_position_vector.size()));

  Config config_out(config_.GetBasis(), relative_position_matrix, atom_vector);

  if (output_folder.empty()) {
    throw std::invalid_argument("Output folder name is empty");
  } else {
    boost::filesystem::create_directories(output_folder);
  }
  config_out.WriteXyzExtended(output_folder + "/" + output_name, auxiliary_lists_, global_info_map);

  return clusters_info_array;
}

std::unordered_set<size_t> SoluteCluster::FindSoluteAtomIndexes() const {
  std::unordered_set<size_t> solute_atoms_hashset;
  for (size_t atom_id = 0; atom_id < config_.GetNumAtoms(); ++atom_id) {
    if (config_.GetElementOfAtom(atom_id) == solvent_element_) {
      solute_atoms_hashset.insert(atom_id);
    }
  }
  return solute_atoms_hashset;
}

std::vector<std::vector<size_t> > SoluteCluster::FindAtomListOfClustersBFSHelper(
    std::unordered_set<size_t> unvisited_atoms_id_set) const {
  std::vector<std::vector<size_t> > cluster_atom_list;
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
      const auto &first_neighbors_list = config_.GetNeighborLists()[0][atom_id];
      for (auto neighbor_id : first_neighbors_list) {
        it = unvisited_atoms_id_set.find(neighbor_id);
        if (it != unvisited_atoms_id_set.end()) {
          visit_id_queue.push(*it);
          unvisited_atoms_id_set.erase(it);
        }
      }
    }
    cluster_atom_list.push_back(atom_list_of_one_cluster);
  }
  return cluster_atom_list;
}

std::vector<std::vector<size_t> > SoluteCluster::FindAtomListOfClusters() const {
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
  // add adjacent solvent neighbors
  for (auto &cluster : cluster_atom_list) {
    std::unordered_set<size_t> cluster_set(cluster.begin(), cluster.end());
    std::unordered_set<size_t> neighbor_set;
    for (auto atom_id : cluster) {
      neighbor_set.insert(atom_id);
      for (auto neighbor_id : config_.GetNeighborLists()[0][atom_id]) {
        neighbor_set.insert(neighbor_id);
      }
    }

    for (auto atom_id : neighbor_set) {
      size_t neighbor_bond_count = 0;
      for (auto neighbor_id : config_.GetNeighborLists()[0][atom_id]) {
        if (cluster_set.find(neighbor_id) != cluster_set.end()
            && config_.GetElementOfAtom(neighbor_id) != solvent_element_) {
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
  std::sort(cluster_atom_list.begin(), cluster_atom_list.end(),
            [](const std::vector<size_t> &a, const std::vector<size_t> &b) {
              return a.size() > b.size();
            });

  return cluster_atom_list;
}
void SoluteCluster::AppendAtomAndLatticeVector(const std::vector<size_t> &cluster_atom_id_list,
                                               std::vector<Element> &atom_vector,
                                               std::vector<Eigen::Vector3d> &relative_position_vector) const {
  for (size_t atom_id : cluster_atom_id_list) {
    atom_vector.push_back(config_.GetAtomVector()[atom_id]);
    relative_position_vector.push_back(config_.GetRelativePositionOfAtom(atom_id));
  }
}

std::map<std::string, size_t> SoluteCluster::GetElementsNumber(
    const std::vector<size_t> &cluster_atom_id_list) const {
  // initialize map with all the element, because some cluster may not have all types of element
  std::map<std::string, size_t> num_atom_in_one_cluster{{"X", 0}};
  for (const auto &element : element_set_) {
    num_atom_in_one_cluster[element.GetElementString()] = 0;
  }
  for (const auto &atom_id : cluster_atom_id_list) {
    num_atom_in_one_cluster.at(config_.GetAtomVector()[atom_id].GetElementString())++;
  }
  return num_atom_in_one_cluster;
}
double SoluteCluster::GetMass(const std::vector<size_t> &cluster_atom_id_list) const {
  double sum_mass = 0;
  for (const auto &atom_id : cluster_atom_id_list) {
    sum_mass += config_.GetAtomVector()[atom_id].GetMass();
  }
  return sum_mass;
}
double SoluteCluster::GetFormationEnergy(const std::vector<size_t> &cluster_atom_id_list) const {
  Config solute_config(solvent_config_);
  double energy_change_solution_to_pure_solvent = 0;
  for (size_t atom_id : cluster_atom_id_list) {
    Element element = config_.GetElementOfAtom(atom_id);
    solute_config.SetElementOfAtom(atom_id, element);
    energy_change_solution_to_pure_solvent += chemical_potential_map_.at(element);
  }
  double energy_change_cluster_to_pure_solvent =
      energy_estimator_.GetEnergyOfCluster(solute_config, cluster_atom_id_list) -
          energy_estimator_.GetEnergyOfCluster(solvent_config_, cluster_atom_id_list);
  return (energy_change_cluster_to_pure_solvent - energy_change_solution_to_pure_solvent);
}
Eigen::Vector3d SoluteCluster::GetGeometryCenter(const std::vector<size_t> &cluster_atom_id_list) const {
  Eigen::Vector3d geometry_center{};
  Eigen::Vector3d sum_cos_theta{};
  Eigen::Vector3d sum_sin_theta{};
  for (size_t atom_id : cluster_atom_id_list) {
    const auto relative_position = config_.GetRelativePositionOfAtom(atom_id);
    for (const int kDim : std::vector<int>{0, 1, 2}) {
      auto theta = relative_position[kDim] * 2 * M_PI;
      sum_cos_theta[kDim] += std::cos(theta);
      sum_sin_theta[kDim] += std::sin(theta);
    }
  }
  auto cos_theta_bar = sum_cos_theta / static_cast<double>(cluster_atom_id_list.size());
  auto sin_theta_bar = sum_sin_theta / static_cast<double>(cluster_atom_id_list.size());
  for (const int kDim : std::vector<int>{0, 1, 2}) {
    double theta_bar = std::atan2(-sin_theta_bar[kDim], -cos_theta_bar[kDim]) + M_PI;
    geometry_center[kDim] = theta_bar / (2 * M_PI);
  }
  return config_.GetBasis() * geometry_center; // Cartesian position
}
Eigen::Vector3d SoluteCluster::GetMassCenter(const std::vector<size_t> &cluster_atom_id_list) const {
  Eigen::Vector3d mass_center{};
  Eigen::Vector3d sum_cos_theta{};
  Eigen::Vector3d sum_sin_theta{};
  double sum_mass = 0;
  for (size_t atom_id : cluster_atom_id_list) {
    const auto relative_position = config_.GetRelativePositionOfAtom(atom_id);
    auto mass = config_.GetAtomVector()[atom_id].GetMass();
    sum_mass += mass;
    for (const int kDim : std::vector<int>{0, 1, 2}) {
      auto theta = relative_position[kDim] * 2 * M_PI;
      sum_cos_theta[kDim] += std::cos(theta) * mass;
      sum_sin_theta[kDim] += std::sin(theta) * mass;
    }
  }
  auto cos_theta_bar = sum_cos_theta / sum_mass;
  auto sin_theta_bar = sum_sin_theta / sum_mass;
  for (const int kDim : std::vector<int>{0, 1, 2}) {
    double theta_bar = std::atan2(-sin_theta_bar[kDim], -cos_theta_bar[kDim]) + M_PI;
    mass_center[kDim] = theta_bar / (2 * M_PI);
  }
  return config_.GetBasis() * mass_center; // Cartesian position
}
Eigen::Matrix3d SoluteCluster::GetMassGyrationTensor(const std::vector<size_t> &cluster_atom_id_list,
                                                     const Eigen::Vector3d &mass_center) const {
  auto relative_mass_center = config_.GetBasis().inverse() * mass_center;
  Eigen::Matrix3d gyration_tensor{};
  double sum_mass = 0;
  for (size_t atom_id : cluster_atom_id_list) {
    const auto relative_position = config_.GetRelativePositionOfAtom(atom_id);
    auto mass = config_.GetAtomVector()[atom_id].GetMass();
    sum_mass += mass;
    for (const int kDim1 : std::vector<int>{0, 1, 2}) {
      auto r1 = relative_position[kDim1] - relative_mass_center[kDim1];
      while (r1 >= 0.5) { r1 -= 1; }
      while (r1 < -0.5) { r1 += 1; }
      for (const int kDim2 : std::vector<int>{0, 1, 2}) {
        auto r2 = relative_position[kDim2] - relative_mass_center[kDim2];
        while (r2 >= 0.5) { r2 -= 1; }
        while (r2 < -0.5) { r2 += 1; }
        gyration_tensor(kDim1, kDim2) += r1 * r2 * mass;
      }
    }
  }
  gyration_tensor /= sum_mass;
  return config_.GetBasis() * config_.GetBasis() * gyration_tensor;
}
Eigen::Matrix3d SoluteCluster::GetMassInertiaTensor(const std::vector<size_t> &cluster_atom_id_list,
                                                    const Eigen::Vector3d &mass_center) const {
  auto relative_mass_center = config_.GetBasis().inverse() * mass_center;
  Eigen::Matrix3d inertia_tensor{};
  for (size_t atom_id : cluster_atom_id_list) {
    auto mass = config_.GetAtomVector()[atom_id].GetMass();
    Eigen::Vector3d relative_distance = config_.GetRelativePositionOfAtom(atom_id) - relative_mass_center;
    for (const int kDim : std::vector<int>{0, 1, 2}) {
      while (relative_distance[kDim] >= 0.5) { relative_distance[kDim] -= 1; }
      while (relative_distance[kDim] < -0.5) { relative_distance[kDim] += 1; }
    }
    auto cartesian_distance = config_.GetBasis() * relative_distance;
    inertia_tensor(0, 0) += mass * (cartesian_distance[1] * cartesian_distance[1]
        + cartesian_distance[2] * cartesian_distance[2]);
    inertia_tensor(0, 1) += -mass * (cartesian_distance[0] * cartesian_distance[1]);
    inertia_tensor(0, 2) += -mass * (cartesian_distance[0] * cartesian_distance[2]);
    inertia_tensor(1, 0) += -mass * (cartesian_distance[1] * cartesian_distance[0]);
    inertia_tensor(1, 1) += mass * (cartesian_distance[2] * cartesian_distance[2]
        + cartesian_distance[0] * cartesian_distance[0]);
    inertia_tensor(1, 2) += -mass * (cartesian_distance[1] * cartesian_distance[2]);
    inertia_tensor(2, 0) += -mass * (cartesian_distance[2] * cartesian_distance[0]);
    inertia_tensor(2, 1) += -mass * (cartesian_distance[2] * cartesian_distance[1]);
    inertia_tensor(2, 2) += mass * (cartesian_distance[0] * cartesian_distance[0]
        + cartesian_distance[1] * cartesian_distance[1]);
  }
  return inertia_tensor;
} // kn

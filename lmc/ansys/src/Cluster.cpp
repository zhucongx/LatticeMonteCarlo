#include "Cluster.h"

#include <queue>
#include <unordered_map>
#include <utility>
#include <filesystem>
#include <cmath>
#include <omp.h>

using json = nlohmann::json;
namespace ansys {
Cluster::Cluster(const cfg::Config &config,
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
    solvent_config_.ChangeAtomElementTypeAtAtom(atom_id, solvent_atom_type);
  }
}

json Cluster::GetClustersInfoAndOutput(
    const std::string &output_folder, const std::string &output_name) {
  auto cluster_to_atom_vector = FindAtomListOfClusters();

  json clusters_info_array = json::array();
  std::vector<cfg::Lattice> lattice_vector;
  std::vector<cfg::Atom> atom_vector;
  for (auto &cluster_atom_id_list: cluster_to_atom_vector) {
    AppendAtomAndLatticeVector(cluster_atom_id_list, atom_vector, lattice_vector);
    json cluster_info = json::object();
    cluster_info["size"] = cluster_atom_id_list.size();
    cluster_info["mass"] = GetMassOfCluster(cluster_atom_id_list);
    cluster_info["geometry_center"] = GetGeometryCenterOfCluster(cluster_atom_id_list);
    cluster_info["mass_center"] = GetMassCenterOfCluster(cluster_atom_id_list);
    cluster_info["elements"] = GetElementNumOfCluster(cluster_atom_id_list);
    cluster_info["energy"] = GetEnergyOfCluster(cluster_atom_id_list);

    clusters_info_array.push_back(cluster_info);
  }
  cfg::Config config_out(config_.GetBasis(), lattice_vector, atom_vector, false);

  if (output_folder.empty()) {
    config_out.WriteConfig(output_name, false);
  } else {
    std::filesystem::create_directories(output_folder);
    config_out.WriteConfig(output_folder + "/" + output_name, false);
  }
  return clusters_info_array;
}

std::unordered_set<size_t> Cluster::FindSoluteAtomIndexes() const {
  std::unordered_set<size_t> solute_atoms_hashset;
  for (const auto &atom: config_.GetAtomVector()) {
    if (atom.GetElement() == solvent_element_) { continue; }
    solute_atoms_hashset.insert(atom.GetId());
  }
  return solute_atoms_hashset;
}

std::vector<std::vector<size_t> > Cluster::FindAtomListOfClustersBFSHelper(
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

std::vector<std::vector<size_t> > Cluster::FindAtomListOfClusters() const {
  auto cluster_atom_list = FindAtomListOfClustersBFSHelper(FindSoluteAtomIndexes());

  // remove small clusters
  auto it = cluster_atom_list.begin();
  while (it != cluster_atom_list.end()) {
    if (it->size() <= smallest_cluster_criteria_) {
      it = cluster_atom_list.erase(it);
    } else {
      ++it;
    }
  }
  std::unordered_set<size_t> all_found_solute_set;
  for (const auto &singe_cluster_vector: cluster_atom_list) {
    std::copy(singe_cluster_vector.begin(),
              singe_cluster_vector.end(),
              std::inserter(all_found_solute_set, all_found_solute_set.end()));
  }
  // add solvent neighbors
  for (const auto &atom: config_.GetAtomVector()) {
    if (atom.GetElement() != solvent_element_)
      continue;
    size_t neighbor_bond_count = 0;
    for (auto neighbor_id: config_.GetFirstNeighborsAtomIdVectorOfAtom(atom.GetId())) {
      const auto &neighbor_type = config_.GetAtomVector()[neighbor_id].GetElement();
      if (neighbor_type != solvent_element_
          && all_found_solute_set.find(neighbor_id) != all_found_solute_set.end())
        neighbor_bond_count++;
    }
    if (neighbor_bond_count >= solvent_bond_criteria_)
      all_found_solute_set.insert(atom.GetId());
  }

  cluster_atom_list = FindAtomListOfClustersBFSHelper(all_found_solute_set);
  return cluster_atom_list;
}
void Cluster::AppendAtomAndLatticeVector(const std::vector<size_t> &cluster_atom_id_list,
                                         std::vector<cfg::Atom> &atom_vector,
                                         std::vector<cfg::Lattice> &lattice_vector) const {
  for (size_t atom_id: cluster_atom_id_list) {
    atom_vector.emplace_back(atom_vector.size(),
                             config_.GetAtomVector()[atom_id].GetElement());
    auto relative_position =
        config_.GetLatticeVector()[config_.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    lattice_vector.emplace_back(lattice_vector.size(),
                                relative_position * config_.GetBasis(), relative_position);
  }
}
std::map<std::string, size_t> Cluster::GetElementNumOfCluster(
    const std::vector<size_t> &cluster_atom_id_list) const {
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
double Cluster::GetMassOfCluster(const std::vector<size_t> &cluster_atom_id_list) const {
  double sum_mass = 0;
  for (const auto &atom_id: cluster_atom_id_list) {
    sum_mass += config_.GetAtomVector()[atom_id].GetElement().GetMass();
  }
  return sum_mass;
}
double Cluster::GetEnergyOfCluster(const std::vector<size_t> &cluster_atom_id_list) const {
  cfg::Config solute_config(solvent_config_);
  double energy_change_solution_to_pure_solvent = 0;
  for (size_t atom_id: cluster_atom_id_list) {
    Element element = config_.GetElementAtAtomId(atom_id);
    solute_config.ChangeAtomElementTypeAtAtom(atom_id, element);
    energy_change_solution_to_pure_solvent += chemical_potential_map_.at(element);
  }
  double energy_change_cluster_to_pure_solvent =
      energy_estimator_.GetEnergyOfCluster(solute_config, cluster_atom_id_list) -
          energy_estimator_.GetEnergyOfCluster(solvent_config_, cluster_atom_id_list);
  return energy_change_cluster_to_pure_solvent - energy_change_solution_to_pure_solvent;
}
Vector_t Cluster::GetGeometryCenterOfCluster(const std::vector<size_t> &cluster_atom_id_list) const {
  Vector_t geometry_center{0, 0, 0};
  Vector_t sum_cos_theta{0, 0, 0};
  Vector_t sum_sin_theta{0, 0, 0};
  for (size_t atom_id: cluster_atom_id_list) {
    auto relative_position =
        config_.GetLatticeVector()[config_.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    for (const auto kDim: All_Dimensions) {
      auto theta = relative_position[kDim] * 2 * M_PI;
      sum_cos_theta[kDim] += std::cos(theta);
      sum_sin_theta[kDim] += std::sin(theta);
    }
  }
  sum_cos_theta /= static_cast<double>(cluster_atom_id_list.size());
  sum_sin_theta /= static_cast<double>(cluster_atom_id_list.size());
  for (const auto kDim: All_Dimensions) {
    double theta_bar = std::atan2(-sum_sin_theta[kDim], -sum_cos_theta[kDim]) + M_PI;
    geometry_center[kDim] = theta_bar / (2 * M_PI);
  }
  return geometry_center;
}
Vector_t Cluster::GetMassCenterOfCluster(const std::vector<size_t> &cluster_atom_id_list) const {
  Vector_t mass_center{0, 0, 0};
  Vector_t sum_cos_theta{0, 0, 0};
  Vector_t sum_sin_theta{0, 0, 0};
  double sum_mass = 0;
  for (size_t atom_id: cluster_atom_id_list) {
    auto relative_position =
        config_.GetLatticeVector()[config_.GetLatticeIdFromAtomId(atom_id)].GetRelativePosition();
    auto mass = config_.GetAtomVector()[atom_id].GetElement().GetMass();
    sum_mass += mass;
    for (const auto kDim: All_Dimensions) {
      auto theta = relative_position[kDim] * 2 * M_PI;
      sum_cos_theta[kDim] += std::cos(theta) * mass;
      sum_sin_theta[kDim] += std::sin(theta) * mass;
    }
  }
  sum_cos_theta /= sum_mass;
  sum_sin_theta /= sum_mass;
  for (const auto kDim: All_Dimensions) {
    double theta_bar = std::atan2(-sum_sin_theta[kDim], -sum_cos_theta[kDim]) + M_PI;
    mass_center[kDim] = theta_bar / (2 * M_PI);
  }
  return mass_center;
}
} // kn

#include "ClustersFinder.h"

#include <queue>
#include <unordered_map>
#include <utility>
#include <filesystem>

namespace ansys {
ClustersFinder::ClustersFinder(std::string cfg_filename,
                               Element solvent_atom_type,
                               size_t smallest_cluster_criteria,
                               size_t solvent_bond_criteria)
    : cfg_filename_(std::move(cfg_filename)),
      config_(cfg::Config::ReadCfg(cfg_filename_)),
      solvent_element_(solvent_atom_type),
      element_set_(config_.GetElementSetWithoutVacancy()),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      solvent_bond_criteria_(solvent_bond_criteria) {}

ClustersFinder::ClusterElementNumMap ClustersFinder::FindClustersAndOutput() {
  // auto cluster_to_atom_map = FindAtomListOfClusters();
  //
  // std::vector<std::map<std::string, size_t> > num_atom_in_clusters_set;
  // for (auto &atom_list: cluster_to_atom_map) {
  //   // initialize map with all the element, because some cluster may not have all types of element
  //   std::map<std::string, size_t> num_atom_in_one_cluster;
  //   for (const auto &element: element_set_) {
  //     num_atom_in_one_cluster[element] = 0;
  //   }
  //
  //   for (const auto &atom_index: atom_list) {
  //     num_atom_in_one_cluster[config_.GetAtomList()[atom_index].GetType()]++;
  //     config_out.AppendAtomWithoutChangingAtomID(config_.GetAtomList()[atom_index]);
  //   }
  //
  //   num_atom_in_clusters_set.push_back(std::move(num_atom_in_one_cluster));
  // }
  // cfg::Config config_out(config_.GetBasis(), config_.GetNumAtoms());
  //
  // std::filesystem::create_directories("cluster");
  // auto output_name("cluster/" + cfg_filename_);
  // auto const pos = output_name.find_last_of('.');
  // output_name.insert(pos, "_cluster");
  // config_out.WriteCfg(output_name, false);
  // return num_atom_in_clusters_set;
}


std::unordered_set<size_t> ClustersFinder::FindSoluteAtomIndexes() const {
  std::unordered_set<size_t> solute_atoms_hashset;
  for (const auto &atom: config_.GetAtomVector()) {
    if (atom.GetElement() == solvent_element_ || atom.GetElement() == ElementName::X) { continue; }
    solute_atoms_hashset.insert(atom.GetId());
  }
  return solute_atoms_hashset;
}

std::vector<std::vector<size_t> > ClustersFinder::FindAtomListOfClustersBFSHelper(
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

      for (auto neighbor_id: config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id)) {
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

std::vector<std::vector<size_t> > ClustersFinder::FindAtomListOfClusters() const {
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

  /*
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
      if (neighbor_type != solvent_element_ && neighbor_type != ElementName::X
          && all_found_solute_set.find(neighbor_id) != all_found_solute_set.end())
        neighbor_bond_count++;
    }
    if (neighbor_bond_count >= solvent_bond_criteria_)
      all_found_solute_set.insert(atom.GetId());
  }

  cluster_atom_list = FindAtomListOfClustersBFSHelper(all_found_solute_set);
  */
  return cluster_atom_list;
}

} // namespace kn

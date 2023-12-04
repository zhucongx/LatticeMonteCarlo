#ifndef LMC_LMC_ANSYS_INCLUDE_CLUSTER_H_
#define LMC_LMC_ANSYS_INCLUDE_CLUSTER_H_
#include <set>
#include <unordered_set>
#include <vector>

#include "Config.h"
#include "EnergyPredictor.h"
#include <nlohmann/json.hpp>

namespace ansys {
class SoluteCluster {
 public:
  SoluteCluster(const cfg::Config &config, Element solvent_atom_type, std::set<Element> element_set,
          size_t smallest_cluster_criteria, size_t solvent_bond_criteria, const pred::EnergyPredictor &energy_estimator,
          const std::map<Element, double> &chemical_potential_map);

  nlohmann::json GetClustersInfoAndOutput(const std::string &output_folder, const std::string &output_name);

 private:
  [[nodiscard]] std::unordered_set<size_t> FindSoluteAtomIndexes() const;
  [[nodiscard]] std::vector<std::vector<size_t>>
  FindAtomListOfClustersBFSHelper(std::unordered_set<size_t> unvisited_atoms_id_set) const;
  // remove smaller clusters and add adjacent atoms
  [[nodiscard]] std::vector<std::vector<size_t>> FindAtomListOfClusters() const;
  void AppendInfoToAuxiliaryListsRepeat(const std::string &key, double value, size_t repeat);
  void AppendAtomAndLatticeVector(const std::vector<size_t> &atom_id_list, std::vector<cfg::Atom> &atom_vector,
                                  std::vector<cfg::Lattice> &lattice_vector) const;
  [[nodiscard]] std::map<std::string, size_t> GetElementsNumber(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] double GetMass(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] double GetFormationEnergy(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] Vector_t GetGeometryCenter(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] Vector_t GetMassCenter(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] Matrix_t GetMassGyrationTensor(const std::vector<size_t> &cluster_atom_id_list,
                                               const Vector_t &mass_center) const;
  [[nodiscard]] Matrix_t GetMassInertiaTensor(const std::vector<size_t> &cluster_atom_id_list,
                                              const Vector_t &mass_center) const;
  const cfg::Config &config_;
  cfg::Config solvent_config_;
  const Element solvent_element_;
  const std::set<Element> element_set_;
  const size_t smallest_cluster_criteria_;
  const size_t solvent_bond_criteria_;
  const pred::EnergyPredictor &energy_estimator_;
  const std::map<Element, double> chemical_potential_map_;
  std::map<std::string, std::vector<double>> auxiliary_lists_{};
};
}    // namespace ansys

#endif    //LMC_LMC_ANSYS_INCLUDE_CLUSTER_H_

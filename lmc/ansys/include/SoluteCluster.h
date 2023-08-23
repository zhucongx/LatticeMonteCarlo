/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 6/14/20 1:27 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/22/23 11:08 PM                                                          *
 **************************************************************************************************/

#ifndef LMC_ANSYS_INCLUDE_SOLUTECLUSTER_H_
#define LMC_ANSYS_INCLUDE_SOLUTECLUSTER_H_
#include <vector>
#include <set>
#include <unordered_set>

#include <nlohmann/json.hpp>
#include "Config.h"
#include "EnergyPredictor.h"

class SoluteCluster {
 public:
  SoluteCluster(const Config &config,
                Element solvent_atom_type,
                std::set<Element> element_set,
                size_t smallest_cluster_criteria,
                size_t solvent_bond_criteria,
                const pred::EnergyPredictor &energy_estimator,
                const std::map<Element, double> &chemical_potential_map);

  nlohmann::json GetClustersInfoAndOutput(const std::string &output_folder,
                                          const std::string &output_name,
                                          const std::map<std::string, Config::ValueVariant> & global_info_map)

 private:
  /*! \brief Query for an unordered set of atom indexes of the solute atoms in the configuration.
   *  @return : An unordered set of atom indexes of the solute atoms.
   */
  [[nodiscard]] std::unordered_set<size_t> FindSoluteAtomIndexes() const;

  /*! \brief A helper function using breadth-first search to find connected clusters of
   *         solute atoms by first nearest neighbor distance.
   *  @return : A 2D vector show all connected atom. Each vector (of atom indexes) represents a cluster.
   */
  [[nodiscard]] std::vector<std::vector<size_t> > FindAtomListOfClustersBFSHelper(
      std::unordered_set<size_t> unvisited_atoms_id_set) const;

  /*! \brief A function to find all clusters.
   *  @return : A 2D vector show all connected atom. Each vector (of atom indexes) represents a cluster.
   */
  [[nodiscard]] std::vector<std::vector<size_t> > FindAtomListOfClusters() const;

  void AppendAtomAndLatticeVector(const std::vector<size_t> &atom_id_list,
                                  std::vector<Element> &atom_vector,
                                  std::vector<Eigen::Vector3d> &relative_position_vector) const;
  [[nodiscard]] std::map<std::string, size_t> GetElementsNumber(
      const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] double GetMass(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] double GetFormationEnergy(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] Eigen::Vector3d GetGeometryCenter(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] Eigen::Vector3d GetMassCenter(const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] Eigen::Matrix3d GetMassGyrationTensor(const std::vector<size_t> &cluster_atom_id_list,
                                                      const Eigen::Vector3d &mass_center) const;
  [[nodiscard]] Eigen::Matrix3d GetMassInertiaTensor(const std::vector<size_t> &cluster_atom_id_list,
                                                     const Eigen::Vector3d &mass_center) const;
  const Config &config_;
  Config solvent_config_;
  const Element solvent_element_;
  const std::set<Element> element_set_;
  const size_t smallest_cluster_criteria_;
  const size_t solvent_bond_criteria_;
  const pred::EnergyPredictor &energy_estimator_;
  const std::map<Element, double> chemical_potential_map_;
  std::map<std::string, Config::VectorVariant> auxiliary_lists_{};
  };

#endif //LMC_ANSYS_INCLUDE_SOLUTECLUSTER_H_

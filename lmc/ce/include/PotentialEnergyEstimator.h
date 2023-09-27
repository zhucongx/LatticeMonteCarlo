/**************************************************************************************************
 * Copyright (c) 2022-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 11/15/22 12:36 PM                                                                       *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/27/23 1:09 PM                                                           *
 **************************************************************************************************/

/*! \file  PotentialEnergyEstimator.h
 *  \brief File for the PotentialEnergyEstimator class definition.
 */
#ifndef LMC_CE_INCLUDE_POTENTIALENERGYESTIMATOR_H_
#define LMC_CE_INCLUDE_POTENTIALENERGYESTIMATOR_H_

#include <set>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>
#include "ClusterExpansion.h"

/*! \brief Class for defining cluster expansion Hamiltonian.
 */
class PotentialEnergyEstimator {
 public:
  PotentialEnergyEstimator(const std::string &predictor_filename,
                           const Config &reference_config,
                           const std::set<Element> &element_set,
                           size_t max_cluster_size,
                           size_t max_bond_order);
  ~PotentialEnergyEstimator();
  /*! \brief Get the encode vector of the configuration, which is the number of appearance of
   *         each cluster types plus void cluster.
   *  \param config : The configuration the code works on
   *  \return       : The encode vector
   */
  [[nodiscard]] std::vector<double> GetEncodeVector(const Config &config) const;
  [[nodiscard]] double GetEnergy(const Config &config) const;
  [[nodiscard]] std::vector<double> GetEncodeVectorOfCluster(
      const Config &config, const std::vector<size_t> &atom_id_list) const;
  [[nodiscard]] double GetEnergyOfCluster(const Config &config,
                                          const std::vector<size_t> &atom_id_list) const;
  [[nodiscard]] std::map<Element, double> GetChemicalPotential(Element solvent_element) const;

 private:
  // const Config &config_;
  const Eigen::VectorXd effective_cluster_interaction_{};
  const std::set<Element> element_set_{};
  const std::set<ClusterType> initialized_cluster_type_set_{};
  const std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> all_lattice_hashset_{};
  const size_t max_cluster_size_{};
  const size_t max_bond_order_{};
};

#endif //LMC_CE_INCLUDE_POTENTIALENERGYESTIMATOR_H_

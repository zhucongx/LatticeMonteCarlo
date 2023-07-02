/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 12/6/21 8:55 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/1/23 12:34 AM                                                           *
 **************************************************************************************************/

/*! \file  Cluster.h
 *  \brief File for the Cluster class definition.
 */

#ifndef LMC_CFG_INCLUDE_CLUSTER_HPP_
#define LMC_CFG_INCLUDE_CLUSTER_HPP_

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <type_traits>
#include <boost/functional/hash.hpp>
#include "Config.h"

/*! \brief Class for defining a cluster of atoms and their positions.
 */
class Cluster {
 public:
  /*! \brief Constructor for setting up the cluster of atoms and their positions.
   *  \param cartesian_position_matrix : The relative position matrix (3，n) of the configuration.
   *  \param atom_vector              : The atom vector of the configuration，n atoms in total.
   */
  Cluster(Eigen::Matrix3Xd cartesian_position_matrix, std::vector<Element> atom_vector);

 private:
  Eigen::Matrix3Xd cartesian_position_matrix_;
  std::vector<Element> element_vector_;
};

#endif //LMC_CFG_INCLUDE_CLUSTER_HPP_

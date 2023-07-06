/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 4:05 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_ANSYS_INCLUDE_SHORTRANGEORDER_H_
#define LMC_ANSYS_INCLUDE_SHORTRANGEORDER_H_

#include "Config.h"

namespace ansys {

class ShortRangeOrder {
 public:
  ShortRangeOrder(const cfg::Config &config, const std::set<Element> &element_set);
  // [[nodiscard]] std::map<std::string, double> FindPairCorrelationCluster(
  //     size_t shell_number, const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] std::map<std::string, double> FindWarrenCowley(size_t shell_number) const;
  [[nodiscard]] std::map<std::string, double> FindProbabilityCluster(
      size_t shell_number, const std::vector<size_t> &cluster_atom_id_list) const;
  [[nodiscard]] std::map<std::string, double> FindProbability(size_t shell_number) const;
 protected:
  const cfg::Config &config_;
  const std::set<Element> &element_set_;
};

} // ansys

#endif //LMC_ANSYS_INCLUDE_SHORTRANGEORDER_H_

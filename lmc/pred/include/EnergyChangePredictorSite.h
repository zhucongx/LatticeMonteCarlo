/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/2/23 10:43 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORSITE_H_
#define LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORSITE_H_
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "AtomCluster.hpp"
#include "EnergyUtility.h"

namespace pred {

class EnergyChangePredictorSite {
 public:
  EnergyChangePredictorSite(const std::string &predictor_filename,
                            const cfg::Config &reference_config,
                            std::set<Element> element_set);
  virtual ~EnergyChangePredictorSite();
  [[nodiscard]] double GetDeFromAtomIdSite(
      const cfg::Config &config, size_t atom_id, Element new_element) const;
 private:
  [[nodiscard]] double GetDeFromLatticeIdSite(
      const cfg::Config &config, size_t lattice_id, Element new_element) const;
  const std::set<Element> element_set_;
  std::vector<double> base_theta_{};
  std::unordered_map<cfg::ElementCluster, size_t, boost::hash<cfg::ElementCluster> >
      initialized_cluster_hashmap_{};

  std::unordered_map<size_t, std::vector<std::vector<std::vector<size_t> > > >
      site_neighbors_hashmap_{};
};

} // pred


#endif //LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORSITE_H_

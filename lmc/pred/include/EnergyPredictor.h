#ifndef LMC_LMC_ANSYS_INCLUDE_ENERGYESTIMATOR_H_
#define LMC_LMC_ANSYS_INCLUDE_ENERGYESTIMATOR_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"

namespace pred {
class EnergyPredictor {
 public:
  EnergyPredictor(const std::string &predictor_filename,
                  
std::set<Element> element_set);
  virtual ~EnergyPredictor();
  [[nodiscard]] double GetEnergy(const cfg::Config &config) const;
  [[nodiscard]] double GetEnergyOfCluster(const cfg::Config &config, const std::vector<size_t> &atom_id_list) const;
  [[nodiscard]] std::map<Element, double> GetChemicalPotential(Element solvent_element) const;

 protected:
  [[nodiscard]] std::vector<double> GetEncode(const cfg::Config &config) const;
  [[nodiscard]] std::vector<double> GetEncodeOfCluster(const cfg::Config &config,
                                                       const std::vector<size_t> &atom_id_list) const;
  const std::set<Element> element_set_;
  Eigen::VectorXd base_theta_{};
  std::unordered_map<cfg::ElementCluster, size_t,
                     boost::hash<cfg::ElementCluster>> initialized_cluster_hashmap_{};
};
}
    // namespace pred
#endif    //LMC_LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

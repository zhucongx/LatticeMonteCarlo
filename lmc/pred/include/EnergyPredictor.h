#ifndef LMC_LMC_ANSYS_INCLUDE_ENERGYESTIMATOR_H_
#define LMC_LMC_ANSYS_INCLUDE_ENERGYESTIMATOR_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"

namespace pred {

class EnergyPredictor {
  public:
    EnergyPredictor(const std::string &predictor_filename,
                    std::set<Element> element_set);
    ~EnergyPredictor();
    [[nodiscard]] std::vector<double> GetEncode(const cfg::Config &config) const;
    [[nodiscard]] std::vector<double> GetEncodeOfCluster(
        const cfg::Config &config, const std::vector<size_t> &atom_id_list) const;

    [[nodiscard]] double GetEnergy(const cfg::Config &config) const;
    [[nodiscard]] double GetEnergyOfCluster(const cfg::Config &config,
                                            const std::vector<size_t> &atom_id_list) const;
    [[nodiscard]] std::map<Element, double> GetChemicalPotential(Element solvent_element) const;
  private:
    std::vector<double> base_theta_{};
    const std::set<Element> element_set_;
    std::unordered_map<cfg::ElementCluster, size_t,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_{};
    // std::vector<cfg::ElementCluster> sorted_cluster_type_vector;
};

} // pred
#endif //LMC_LMC_ANSYS_INCLUDE_ENERGYESTIMATOR_H_

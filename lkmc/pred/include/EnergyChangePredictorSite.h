#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORSITE_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORSITE_H_
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"

namespace pred {

class EnergyChangePredictorSite {
  public:
    EnergyChangePredictorSite(const std::string &predictor_filename,
                              const cfg::Config &reference_config,
                              std::set<Element> element_set);
    virtual ~EnergyChangePredictorSite();
    [[nodiscard]] double GetDiffFromAtomId(
        const cfg::Config &config, size_t atom_id, Element new_element) const;
  private:
    [[nodiscard]] double GetDiffFromLatticeIdPair(
        const cfg::Config &config, size_t lattice_id, Element new_element) const;
    const std::set<Element> element_set_;
    std::vector<double> base_theta_{};
    std::unordered_map<cfg::ElementCluster, size_t,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_;

    std::unordered_map<size_t, std::vector<std::vector<std::vector<size_t> > > >
        site_neighbors_hashmap_;
};

} // namespace pred


#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORSITE_H_

#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORSTATE_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORSTATE_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"

namespace pred {

class EnergyPredictorState {
  public:
    EnergyPredictorState(const std::string &predictor_filename,
                    const cfg::Config &reference_config,
                    const std::set<Element> &type_set);
    virtual ~EnergyPredictorState();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
  protected:
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    const std::vector<std::vector<std::vector<size_t> > > mapping_state_{};
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_state_hashmap_;
    std::vector<double> base_theta_{};
    std::unordered_map<Element,
                       std::vector<double>,
                       boost::hash<Element> > element_theta_;
    std::unordered_map<Element,
                       std::unordered_map<cfg::ElementCluster, size_t,
                                          boost::hash<cfg::ElementCluster> >,
                       boost::hash<Element> > element_initialized_cluster_hashmap_;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORSTATE_H_

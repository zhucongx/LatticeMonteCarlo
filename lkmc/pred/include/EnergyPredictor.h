#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "ElementCluster.hpp"
namespace pred {

class EnergyPredictor {
  public:
    EnergyPredictor(const std::string &predictor_filename,
                    const cfg::Config &reference_config,
                    const std::set<Element> &type_set);
    virtual ~EnergyPredictor();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
  protected:
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    const std::vector<std::vector<std::vector<size_t> > > cluster_mapping_{};
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_hashmap_;
    std::vector<double> base_theta_{};
    std::unordered_map<Element,
                       std::vector<double>,
                       boost::hash<Element> > element_theta_;
    std::unordered_map<Element,
                       std::unordered_map<cfg::ElementCluster, int,
                                          boost::hash<cfg::ElementCluster> >,
                       boost::hash<Element> > element_initialized_cluster_hashmap_;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

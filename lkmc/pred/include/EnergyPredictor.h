#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
namespace pred {
std::unordered_map<
    cfg::ElementCluster, int, boost::hash<cfg::ElementCluster> > InitializeClusterHashMap(
    const std::set<Element> &type_set);
std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingState(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingState(
    const cfg::Config &config);
std::vector<cfg::Lattice> GetSortedLatticeVectorState(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::array<std::vector<double>, 2> GetEncodesFromMapState(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    const std::unordered_map<cfg::ElementCluster,
                             int,
                             boost::hash<cfg::ElementCluster> > &initialized_cluster_hashmap,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping);

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

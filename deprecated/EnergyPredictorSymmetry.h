#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORSYMMETRY_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORSYMMETRY_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "EnergyUtility.h"
namespace pred {

class EnergyPredictorSymmetry {
  public:
    EnergyPredictorSymmetry(const std::string &predictor_filename,
                            const cfg::Config &reference_config,
                            const std::set<Element> &type_set);
    virtual ~EnergyPredictorSymmetry();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
  protected:
    [[nodiscard]] double GetE0(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair,
                               Element migration_element) const;
    [[nodiscard]] double GetE(const cfg::Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair,
                              Element migration_element) const;
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  private:
    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mm2_;
    std::unordered_map<Element,
                       ParametersE0,
                       boost::hash<Element>> e0_element_parameters_hashmap_;
    std::unordered_map<Element,
                       ParametersDE,
                       boost::hash<Element>> dE_element_parameters_hashmap_;

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mm2_hashmap_;
    std::set<Element> type_set_;

};
} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORSYMMETRY_H_

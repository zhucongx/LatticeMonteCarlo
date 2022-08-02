#ifndef LKMC_LKMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_
#define LKMC_LKMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"

namespace pred {

class VacancyMigrationPredictorQuartic {
  public:
    VacancyMigrationPredictorQuartic(const std::string &predictor_filename,
                                     const cfg::Config &reference_config,
                                     std::set<Element> element_set);
    virtual ~VacancyMigrationPredictorQuartic();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
  protected:
    [[nodiscard]] double GetDe(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetKs(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetD(const cfg::Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

    const std::set<Element> element_set_;
    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mm2_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_state_;

    std::vector<double> base_theta_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mm2_hashmap_;
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_state_hashmap_;

    std::unordered_map<Element,
                       ParametersQuartic,
                       boost::hash<Element>> element_parameters_hashmap_;

    std::unordered_map<cfg::ElementCluster, size_t,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_;

};
} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_

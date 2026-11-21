#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include <Eigen/Dense>
#include "Config.h"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"

namespace pred {
class VacancyMigrationPredictorQuartic {
  public:
    VacancyMigrationPredictorQuartic(const std::string &predictor_filename,
                                     const cfg::Config &reference_config,
                                     std::set<Element> element_set);
    virtual ~VacancyMigrationPredictorQuartic();
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  protected:
    [[nodiscard]] double GetDe(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetKs(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetD(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] size_t GetPairFlatIndex(const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

    const std::set<Element> element_set_;
    const OneHotEncoder encoder_mmm_;
    const OneHotEncoder encoder_mm2_;
    const FlatClusterMapping mapping_state_;
    const size_t num_pair_slots_{};

    std::vector<std::vector<size_t>> site_bond_cluster_state_flat_{};
    std::vector<std::vector<size_t>> site_bond_cluster_mmm_flat_{};
    std::vector<std::vector<size_t>> site_bond_cluster_mm2_flat_{};
    // The mmm and mm2 lists of a pair hold the sites of its state list in another order, which is the same for every
    // pair: entry n of such a list is entry positions[n] of the state list. `mm2_backward_positions_` is that order
    // for the mm2 list of the reversed pair.
    std::vector<size_t> mmm_positions_{};
    std::vector<size_t> mm2_positions_{};
    std::vector<size_t> mm2_backward_positions_{};

    Eigen::VectorXd base_theta_{};
    std::unordered_map<Element, ParametersQuartic, boost::hash<Element>> element_parameters_hashmap_{};
    std::vector<std::unordered_map<size_t, size_t>> neighbor_slot_lookup_;

    std::unordered_map<cfg::ElementCluster, size_t,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_{};
    std::map<cfg::ElementCluster, int> ordered_map_{};
    ClusterIndexer state_cluster_indexer_{};
};
} // pred
#endif //LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_

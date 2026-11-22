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
// The barrier of a jump as the fitted quartic form a x^4 + b x^3 + c x^2 of the energy along its path gives it, from
// the energy change dE of the jump, the distance D between the two minima of the path and the sum Ks of the
// curvatures at them. With u = 16 dE / (Ks D^2) it is (27 / 16384) Ks D^2 (u + 4/3)^3 (4 - u).
//
// The form has two minima only for |dE| < Ks D^2 / 12. At dE = -Ks D^2 / 12 the minimum the jump starts from merges
// with the saddle point: the returned barrier is zero there and negative below. At dE = Ks D^2 / 12 the same happens
// to the minimum the jump ends in: the returned barrier is dE there and smaller than dE above, so that the jump back
// has a negative barrier.
[[nodiscard]] double GetQuarticModelBarrier(double dE, double D, double Ks);
// The barrier the predictor returns: that of the model, but never less than the energy the jump costs and never less
// than zero, max(model, dE, 0). The two limits take effect exactly where the model has no barrier forward or
// backward, so a barrier the model does have is returned unchanged. The barriers of a jump and of the jump back still
// differ by dE, which keeps the rates in detailed balance.
[[nodiscard]] double GetQuarticBarrier(double dE, double D, double Ks);

class VacancyMigrationPredictorQuartic {
  public:
    VacancyMigrationPredictorQuartic(const std::string &predictor_filename,
                                     const cfg::Config &reference_config,
                                     std::set<Element> element_set);
    virtual ~VacancyMigrationPredictorQuartic();
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
    // The barrier is GetQuarticBarrier() of the jump.
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    // GetQuarticModelBarrier() of the jump: the barrier of the model as it is, which tells how far outside the range
    // of the quartic form a jump without a barrier lies. It is evaluated anew on every call.
    [[nodiscard]] double GetModelBarrierFromLatticeIdPair(
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

    // The sites around every pair in the order of GetSortedLatticeVectorStateOfPair(), its state list.
    std::vector<std::vector<size_t>> site_bond_cluster_state_flat_{};
    // The mmm and mm2 lists of a pair hold the sites of its state list in another order, which is the same for every
    // pair: entry n of such a list is entry positions[n] of the state list. `mm2_backward_positions_` is that order
    // for the mm2 list of the reversed pair. The lists themselves are therefore not stored.
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

#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"
#include <Eigen/Dense> // Added for Eigen::VectorXd

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
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  protected:
    [[nodiscard]] double GetDe(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetKs(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetD(const cfg::Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] size_t GetPairFlatIndex(const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    const std::set<Element> element_set_;
    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mm2_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_state_;
    // Total slots in flattened neighbor storage: num_sites * NN per site.
    const size_t num_pair_slots_;
    // Per-site map: neighbor lattice id -> neighbor slot index (0..11).
    // Built once from the reference config to avoid per-call linear scans.
    std::vector<std::unordered_map<size_t, size_t> > neighbor_slot_lookup_;

    Eigen::VectorXd base_theta_{}; // Changed from std::vector<double>

    // Flattened storage of site-bond clusters to remove hot-path unordered_map lookups.
    // Indexing: flat_idx = site * constants::kNumFirstNearestNeighbors + neighbor_slot.
    std::vector<std::vector<size_t> > site_bond_cluster_mmm_flat_{};
    std::vector<std::vector<size_t> > site_bond_cluster_mm2_flat_{};
    std::vector<std::vector<size_t> > site_bond_cluster_state_flat_{};

    std::unordered_map<Element,
                       ParametersQuartic,
                       boost::hash<Element>> element_parameters_hashmap_{};

    std::unordered_map<cfg::ElementCluster, size_t,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_{};
    std::map<cfg::ElementCluster, int> ordered_map_{};
    ClusterIndexer state_cluster_indexer_{};
};
} // pred
#endif //LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTIC_H_

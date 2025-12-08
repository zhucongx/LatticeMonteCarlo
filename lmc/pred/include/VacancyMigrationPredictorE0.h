#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORE0_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORE0_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"
#include <Eigen/Dense>

namespace pred {
class VacancyMigrationPredictorE0 {
  public:
    VacancyMigrationPredictorE0(const std::string &predictor_filename,
                                const cfg::Config &reference_config,
                                std::set<Element> element_set);
    virtual ~VacancyMigrationPredictorE0();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  protected:
    [[nodiscard]] double GetDe(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetE0(const cfg::Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    const std::set<Element> element_set_;
    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_state_;

    Eigen::VectorXd base_theta_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_state_hashmap_{};

    std::unordered_map<Element,
                       ParametersE0,
                       boost::hash<Element>> element_parameters_hashmap_{};

    std::unordered_map<cfg::ElementCluster, size_t,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_{};
    std::map<cfg::ElementCluster, int> ordered_map_{};

};
} // pred
#endif //LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORE0_H_

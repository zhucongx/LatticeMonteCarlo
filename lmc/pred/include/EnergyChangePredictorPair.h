#ifndef LMC_LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORPAIR_H_
#define LMC_LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORPAIR_H_

#include <string>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "EnergyUtility.h"
#include <Eigen/Dense>

namespace pred {
class EnergyChangePredictorPair {
  public:
    EnergyChangePredictorPair(const std::string &predictor_filename,
                             const cfg::Config &reference_config,
                             std::set<Element> element_set);
    virtual ~EnergyChangePredictorPair();
    [[nodiscard]] double GetDeFromAtomIdPair(const cfg::Config &config,
                                             const std::pair<size_t, size_t> &atom_id_jump_pair) const;
    [[nodiscard]] double GetDeFromLatticeIdPair(const cfg::Config &config,
                                                const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  private:
    const std::set<Element> element_set_;
    const std::vector<std::vector<std::vector<size_t>>> bond_mapping_state_;

    Eigen::VectorXd base_theta_{};
    ClusterIndexer cluster_indexer_{};

    std::unordered_map<std::pair<size_t, size_t>, std::vector<size_t>, boost::hash<std::pair<size_t, size_t>>>
        bond_state_hashmap_{};
};
} // pred

#endif //LMC_LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORPAIR_H_

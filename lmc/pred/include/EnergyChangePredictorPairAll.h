#ifndef LMC_LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORPAIRALL_H_
#define LMC_LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORPAIRALL_H_
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyUtility.h"

namespace pred {
class EnergyChangePredictorPairAll {
  public:
    EnergyChangePredictorPairAll(const std::string &predictor_filename,
                                 const cfg::Config &reference_config,
                                 std::set<Element> element_set);
    virtual ~EnergyChangePredictorPairAll();
    [[nodiscard]] double GetDeFromAtomIdPair(
        const cfg::Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair) const;
    [[nodiscard]] double GetDeFromLatticeIdPair(
        const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  private:
    [[nodiscard]] double GetDeFromLatticeIdPairWithCoupling(
        const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetDeFromLatticeIdPairWithoutCoupling(
        const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetDeFromLatticeIdSite(
        const cfg::Config &config, size_t lattice_id, Element new_element) const;
    [[nodiscard]] double GetDeHelper(
        const std::unordered_map<cfg::ElementCluster, size_t,
                                 boost::hash<cfg::ElementCluster> > &start_hashmap,
        const std::unordered_map<cfg::ElementCluster, size_t,
                                 boost::hash<cfg::ElementCluster> > &end_hashmap,
        const std::map<cfg::ElementCluster, int> &ordered) const;

    const std::set<Element> element_set_;
    const std::vector<std::vector<std::vector<size_t> > > site_mapping_state_;

    std::vector<double> base_theta_{};
    std::unordered_map<cfg::ElementCluster, size_t,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_;

    std::unordered_map<size_t,
                       std::vector<size_t> > site_state_hashmap_;
    std::unordered_map<size_t,
                       std::unordered_set<size_t> > neighboring_sites_hashmap_;
};
} // pred
#endif //LMC_LMC_PRED_INCLUDE_ENERGYCHANGEPREDICTORPAIRALL_H_

#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DESTATE_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DESTATE_H_
#include "EnergyPredictor.h"
#include "ElementCluster.hpp"
namespace pred {

class EnergyPredictorE0DEState : public EnergyPredictor {
  public:
    EnergyPredictorE0DEState(const std::string &predictor_filename,
                             const cfg::Config &reference_config,
                             const std::set<Element> &type_set);
    ~EnergyPredictorE0DEState() override;
  protected:
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const override;

  private:
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<std::vector<std::vector<size_t> > >,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_mapping_hashmap_;
    std::vector<double> base_theta_{};
    std::unordered_map<Element,
                       std::vector<double>,
                       boost::hash<Element> > element_theta_;
    std::unordered_map<Element,
                       std::unordered_map<ElementCluster, int,
                                          boost::hash<ElementCluster> >,
                       boost::hash<Element> > element_initialized_cluster_hashmap_;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DESTATE_H_

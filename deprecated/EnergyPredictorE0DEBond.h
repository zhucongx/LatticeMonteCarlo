#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DEBOND_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DEBOND_H_
#include "EnergyPredictor.h"
#include "ElementBond.hpp"

namespace pred {

class EnergyPredictorE0DEBond : public EnergyPredictor {
  public:
    EnergyPredictorE0DEBond(const std::string &predictor_filename,
                            const cfg::Config &reference_config,
                            const std::set<Element> &type_set);
    ~EnergyPredictorE0DEBond() override;
  protected:
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const override;
    [[nodiscard]] double GetE0(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair,
                               Element migration_element) const;
    [[nodiscard]] double GetDE(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  private:
    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;

    std::unordered_map<Element,
                       ParametersE0,
                       boost::hash<Element>> element_parameters_hashmap_;
    std::vector<double> theta_bond_{};
    const std::unordered_map<ElementBond, int,
                             boost::hash<ElementBond> > initialized_bond_hashmap_;
};

} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DEBOND_H_

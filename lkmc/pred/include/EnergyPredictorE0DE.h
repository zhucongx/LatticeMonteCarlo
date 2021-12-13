#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DE_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DE_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "EnergyPredictor.h"
namespace pred {

class EnergyPredictorE0DE : public EnergyPredictor {
  public:
    EnergyPredictorE0DE(const std::string &predictor_filename,
                    const cfg::Config &reference_config,
                    const std::set<Element> &type_set);
    ~EnergyPredictorE0DE() override;
  protected:
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const override;
    [[nodiscard]] double GetE0(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair,
                               Element migration_element) const;
    [[nodiscard]] double GetE(const cfg::Config &config,
                              size_t lattice_id,
                              Element migration_element) const;
  private:
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_{};
    const std::vector<std::vector<std::vector<size_t> > > mapping_periodic_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;
    std::unordered_map<size_t,
                       std::vector<size_t> > site_cluster_periodic_hashmap_;

    std::unordered_map<Element,
                       ParametersE0,
                       boost::hash<Element>> element_parameters_hashmap_;
    ParametersDE parameters_dE_;
};

} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORE0DE_H_

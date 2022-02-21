#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "ClusterExpansionPredictor.h"
namespace pred {
class EnergyPredictor {
  public:
    EnergyPredictor();
    virtual ~EnergyPredictor();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
  protected:
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const = 0;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

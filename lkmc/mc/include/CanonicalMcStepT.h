#ifndef LKMC_LKMC_MC_INCLUDE_CANONICALMCSTEPT_H_
#define LKMC_LKMC_MC_INCLUDE_CANONICALMCSTEPT_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "McAbstract.h"
#include "EnergyChangePredictorPairAll.h"
namespace mc {
class CanonicalMcStepT : public McAbstract {
  public:
    CanonicalMcStepT(cfg::Config config,
                     unsigned long long int log_dump_steps,
                     unsigned long long int config_dump_steps,
                     unsigned long long int maximum_steps,
                     unsigned long long int thermodynamic_averaging_steps,
                     double initial_temperature,
                     double decrement_temperature,
                     const std::set<Element> &element_set,
                     const std::string &json_coefficients_filename);
    void Simulate() override;
  private:
    inline void UpdateTemperature();
    inline void Dump();
    std::pair<size_t, size_t> GenerateLatticeIdJumpPair();

    // simulation parameters
    const double initial_temperature_;
    const double decrement_temperature_;

    // helpful properties
    const pred::EnergyChangePredictorPairAll energy_predictor_;
    mutable std::uniform_int_distribution<size_t> atom_index_selector_;
};
} // mc

#endif //LKMC_LKMC_MC_INCLUDE_CANONICALMCSTEPT_H_

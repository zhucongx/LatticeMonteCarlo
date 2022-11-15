#ifndef LKMC_LKMC_ANSYS_INCLUDE_CANONICALMCNEARSTEPT_H_
#define LKMC_LKMC_ANSYS_INCLUDE_CANONICALMCNEARSTEPT_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "EnergyChangePredictorPair.h"
#include "EnergyChangePredictorSite.h"
namespace ansys {
class CanonicalMcNearStepT {
  public:
    CanonicalMcNearStepT(cfg::Config config,
                         Element solvent_element,
                         const std::set<Element> &solute_element_set,
                         unsigned long long int log_dump_steps,
                         unsigned long long int config_dump_steps,
                         unsigned long long int maximum_number,
                         double initial_temperature,
                         double decrement_temperature,
                         const std::string &json_coefficients_filename);
    void Simulate();
  private:
    inline void UpdateTemperature();
    inline void Dump(std::ofstream &ofs);
    std::pair<size_t, size_t> GenerateAtomIdJumpPair();
    // simulation parameters
    cfg::Config config_;
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_number_;
    // simulation statistics
    unsigned long long int steps_{0};
    double energy_{0.0};
    const double initial_temperature_;
    const double decrement_temperature_;
    double temperature_;
    double beta_;

    const pred::EnergyChangePredictorPair energy_predictor_;
    mutable std::mt19937_64 generator_;
    mutable std::uniform_int_distribution<size_t> atom_index_selector_;
    mutable std::uniform_int_distribution<size_t> neighbor_index_selector_;
    mutable std::uniform_real_distribution<double> one_distribution_{0.0, 1.0};
};
} // namespace ansys

#endif //LKMC_LKMC_ANSYS_INCLUDE_CANONICALMCNEARSTEPT_H_

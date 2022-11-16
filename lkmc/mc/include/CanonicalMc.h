#ifndef LKMC_LKMC_MC_INCLUDE_CANONICALMC_H_
#define LKMC_LKMC_MC_INCLUDE_CANONICALMC_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "EnergyChangePredictorPair.h"
#include "EnergyChangePredictorSite.h"
namespace mc {
class CanonicalMc {
  public:
    CanonicalMc(cfg::Config config,
                const std::set<Element> &element_set,
                unsigned long long int log_dump_steps,
                unsigned long long int config_dump_steps,
                unsigned long long int maximum_number,
                double temperature,
                const std::string &json_coefficients_filename);
    void Simulate();
  private:
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
    const double beta_;

    const pred::EnergyChangePredictorPair energy_predictor_;
    mutable std::mt19937_64 generator_;
    mutable std::uniform_int_distribution<size_t> atom_index_selector_;
    mutable std::uniform_int_distribution<size_t> neighbor_index_selector_;
    mutable std::uniform_real_distribution<double> one_distribution_{0.0, 1.0};
};
} // mc

#endif //LKMC_LKMC_MC_INCLUDE_CANONICALMC_H_

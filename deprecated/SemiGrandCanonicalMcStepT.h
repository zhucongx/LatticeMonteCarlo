#ifndef LKMC_LKMC_MC_INCLUDE_SEMIGRANDCANONICALMCSTEPT_H_
#define LKMC_LKMC_MC_INCLUDE_SEMIGRANDCANONICALMCSTEPT_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "EnergyChangePredictorSite.h"
namespace mc {
class SemiGrandCanonicalMcStepT {
  public:
    SemiGrandCanonicalMcStepT(cfg::Config config,
                              const std::set<Element> &element_set,
                              unsigned long long int log_dump_steps,
                              unsigned long long int config_dump_steps,
                              unsigned long long int maximum_steps,
                              unsigned long long int thermodynamic_averaging_steps,
                              double initial_temperature,
                              double decrement_temperature,
                              const std::string &json_coefficients_filename);
    void Simulate();
  private:
    inline void UpdateTemperature();
    inline void Dump(std::ofstream &ofs);
    std::pair<size_t, Element> GenerateAtomIdChangeSite();
    // simulation parameters
    cfg::Config config_;
    const std::vector<Element> element_vector_;
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_steps_;
    // simulation statistics
    unsigned long long int steps_{0};
    double energy_{0.0};
    const double initial_temperature_;
    const double decrement_temperature_;
    double temperature_;
    double beta_;

    const pred::EnergyChangePredictorSite energy_predictor_;
    std::unordered_map<Element, double, boost::hash<Element>> chemical_potential_;

    mutable std::mt19937_64 generator_;
    mutable std::uniform_int_distribution<size_t> atom_index_selector_;
    mutable std::uniform_int_distribution<size_t> element_index_selector_;
    mutable std::uniform_real_distribution<double> one_distribution_{0.0, 1.0};
};
} // mc

#endif //LKMC_LKMC_MC_INCLUDE_SEMIGRANDCANONICALMCSTEPT_H_

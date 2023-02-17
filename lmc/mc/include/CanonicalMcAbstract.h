#ifndef LMC_LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "McAbstract.h"
#include "EnergyChangePredictorPairAll.h"
namespace mc {

class CanonicalMcAbstract : public McAbstract {
  public:
    CanonicalMcAbstract(cfg::Config config,
                        unsigned long long int log_dump_steps,
                        unsigned long long int config_dump_steps,
                        unsigned long long int maximum_steps,
                        unsigned long long int thermodynamic_averaging_steps,
                        double temperature,
                        // double initial_temperature,
                        // double decrement_temperature,
                        const std::set<Element> &element_set,
                        const std::string &json_coefficients_filename);
    void Simulate() override = 0;
  protected:
    // void UpdateTemperature();
    void Dump();
    std::pair<size_t, size_t> GenerateLatticeIdJumpPair();
    void SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair, double dE);
    // simulation parameters
    // const double initial_temperature_;
    // const double decrement_temperature_;
    // helpful properties
    const pred::EnergyChangePredictorPairAll energy_change_predictor_;
    mutable std::uniform_int_distribution<size_t> atom_index_selector_;
};

} // mc

#endif //LMC_LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_

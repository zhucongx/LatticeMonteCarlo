#ifndef LMC_LMC_MC_INCLUDE_KINETICMCSECONDOMPI_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCSECONDOMPI_H_
#include <omp.h>
#include "KineticMcSecondMpi.h"
namespace mc {
//  j -> k -> i -> l
//       |
// current position

// KineticMcSecondMpi with the twelve jumps i -> l of a process spread over its OpenMP threads.
class KineticMcSecondOmpi : public KineticMcSecondMpi {
  public:
    KineticMcSecondOmpi(cfg::Config config,
                        unsigned long long int log_dump_steps,
                        unsigned long long int config_dump_steps,
                        unsigned long long int maximum_steps,
                        unsigned long long int thermodynamic_averaging_steps,
                        unsigned long long int restart_steps,
                        double restart_energy,
                        double restart_time,
                        double temperature,
                        const std::set<Element> &element_set,
                        const std::string &json_coefficients_filename,
                        const std::string &time_temperature_filename,
                        bool is_early_stop,
                        bool is_solute_disp,
                        unsigned long long int seed);
  protected:
    void EvaluateSecondJumps(size_t k_lattice_id, size_t i_lattice_id) override;
};

} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCSECONDOMPI_H_

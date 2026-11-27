#ifndef LMC_LMC_MC_INCLUDE_KINETICMCSECONDSERIAL_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCSECONDSERIAL_H_
#include <array>
#include <omp.h>
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
namespace mc {
//  j -> k -> i -> l
//       |
// current position

// Second-order KMC in one process: the twelve first jumps k -> i, and the twelve jumps i -> l that can follow each of
// them, are evaluated one after the other.
class KineticMcSecondSerial : public KineticMcSecondAbstract {
  public:
    KineticMcSecondSerial(cfg::Config config,
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
    ~KineticMcSecondSerial() override;
  protected:
    void BuildEventList() override;
    double CalculateTime() override;
    // For every first jump k -> i, the sum of the rates of the twelve jumps i -> l that can follow it.
    std::array<double, kEventListSize> total_rate_i_list_{};
};

} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCSECONDSERIAL_H_

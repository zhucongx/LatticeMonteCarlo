#ifndef LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "VacancyMigrationPredictorQuarticLru.h"
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
namespace mc {
//  j -> k -> i -> l
//       |
// current position

class KineticMcChainOmpi : public KineticMcChainAbstract {
  public:
    KineticMcChainOmpi(cfg::Config config,
                       unsigned long long int log_dump_steps,
                       unsigned long long int config_dump_steps,
                       unsigned long long int maximum_steps,
                       unsigned long long int thermodynamic_averaging_steps,
                       unsigned long long int restart_steps,
                       double restart_energy,
                       double restart_time,
                       double temperature,
                       const std::set<Element> &element_set,
                       const std::string &json_coefficients_filename);
  protected:
    void BuildEventList() override;
    double CalculateTime() override;
};

} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_

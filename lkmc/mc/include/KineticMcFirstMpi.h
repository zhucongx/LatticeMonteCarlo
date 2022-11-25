#ifndef LKMC_LKMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_
#define LKMC_LKMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_
#include <mpi.h>
#include "KineticMcAbstract.h"
namespace mc {
class KineticMcFirstMpi : public KineticMcAbstract {
  public:
    KineticMcFirstMpi(cfg::Config config,
                      unsigned long long int log_dump_steps,
                      unsigned long long int config_dump_steps,
                      unsigned long long int maximum_steps,
                      unsigned long long int thermodynamic_averaging_steps,
                      double temperature,
                      const std::set<Element> &element_set,
                      unsigned long long int restart_steps,
                      double restart_energy,
                      double restart_time,
                      const std::string &json_coefficients_filename);
    ~KineticMcFirstMpi() override;
    void Simulate() override;
  protected:
    void BuildEventList();
    double CalculateTime();
    // helpful properties
    int world_rank_{-1};
};
} // mc

#endif //LKMC_LKMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_

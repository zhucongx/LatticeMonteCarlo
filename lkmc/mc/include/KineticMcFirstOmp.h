#ifndef LKMC_LKMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#define LKMC_LKMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#include <random>
#include <omp.h>
#include "VacancyMigrationPredictorQuarticLru.h"
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
namespace mc {
class KineticMcFirstOmp : public KineticMcAbstract {
  public:
    KineticMcFirstOmp(cfg::Config config,
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
    ~KineticMcFirstOmp() override;
    void Simulate() override;
  protected:
    void BuildEventList();
    double CalculateTime();
};

} // mc

#endif //LKMC_LKMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_

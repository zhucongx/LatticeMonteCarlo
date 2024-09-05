#ifndef LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#include <random>
#include <omp.h>
#include "VacancyMigrationPredictorQuarticLru.h"
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
namespace mc {
class KineticMcFirstOmp : public KineticMcFirstAbstract {
  public:
    KineticMcFirstOmp(cfg::Config config,
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
                      bool is_rate_corrector,
                      const Vector_t &vacancy_trajectory,
                      bool is_early_stop);
    ~KineticMcFirstOmp() override;
  protected:
    void BuildEventList() override;
    double CalculateTime() override;
};
} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_

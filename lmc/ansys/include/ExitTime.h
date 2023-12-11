#ifndef LMC_LMC_ANSYS_INCLUDE_EXITTIME_H
#define LMC_LMC_ANSYS_INCLUDE_EXITTIME_H

#include "Config.h"
#include "VacancyMigrationPredictorQuartic.h"
namespace ansys {
class ExitTime {
 public:
  ExitTime(const cfg::Config &config,
           const Element &solvent_element,
           const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor,
           double temperature);
  [[nodiscard]] std::pair<std::vector<std::vector<double>>, std::vector<double>> GetBarrierListAndExitTime() const;
 protected:
  const cfg::Config &config_;
  const Element &solvent_element_;
  const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor_;
  const double beta_;
};

}    // namespace ansys

#endif    //LMC_LMC_ANSYS_INCLUDE_EXITTIME_H

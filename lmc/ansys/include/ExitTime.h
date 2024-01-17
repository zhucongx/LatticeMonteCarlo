#ifndef LMC_LMC_ANSYS_INCLUDE_EXITTIME_H
#define LMC_LMC_ANSYS_INCLUDE_EXITTIME_H

#include "Config.h"
#include "EnergyChangePredictorPairSite.h"
#include "VacancyMigrationPredictorQuartic.h"

namespace ansys {
class ExitTime {
 public:
  ExitTime(const cfg::Config &config,
           const Element &solvent_element,
           double temperature,
           const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor,
           const pred::EnergyChangePredictorPairSite &energy_change_predictor_site,
           const std::map<Element, double> &chemical_potential);
  [[nodiscard]] std::tuple<std::vector<std::vector<double>>,  std::vector<double>, std::vector<double>> GetBarrierListAndExitTime() const;
  [[nodiscard]] std::vector<double> GetAverageBarriers(const std::unordered_set<size_t> &atom_id_set) const;
  // [[nodiscard]] double GetLocalBindingEnergy() const;
  [[nodiscard]] std::map<Element, std::vector<double>> GetBindingEnergy() const;
  // [[nodiscard]] std::vector<double> GetProfileEnergy() const;

 protected:
  const cfg::Config &config_;
  const Element solvent_element_;
  const double beta_;
  const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor_;
  const pred::EnergyChangePredictorPairSite &energy_change_predictor_pair_site_;
  const std::map<Element, double> chemical_potential_;
};
}    // namespace ansys

#endif    //LMC_LMC_ANSYS_INCLUDE_EXITTIME_H
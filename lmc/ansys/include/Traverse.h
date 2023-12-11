#ifndef LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#define LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#include "EnergyPredictor.h"
#include "ExitTime.h"
#include "ShortRangeOrder.h"
#include "SoluteCluster.h"
#include "VacancyMigrationPredictorQuartic.h"

namespace ansys {

class Traverse {
 public:
  Traverse(unsigned long long int initial_steps,
           unsigned long long int increment_steps,
           Element solvent_element,
           std::set<Element> element_set,
           size_t smallest_cluster_criteria,
           size_t solvent_bond_criteria,
           const std::string &predictor_filename,
           std::string log_type,
           std::string config_type);
  virtual ~Traverse();
  void RunAnsys() const;
  void RunReformat() const;

 private:
  const unsigned long long initial_steps_;
  const unsigned long long increment_steps_;
  unsigned long long final_number_;
  Element solvent_element_;
  const std::set<Element> element_set_;
  size_t smallest_cluster_criteria_;
  size_t solvent_bond_criteria_;
  std::unordered_map<unsigned long long, double> filename_energy_hashset_{}, filename_temperature_hashset_{},
      filename_time_hashset_{};
  const std::string log_type_;
  const std::string config_type_;
  const pred::EnergyPredictor energy_estimator_;
  const pred::VacancyMigrationPredictorQuartic vacancy_migration_predictor_;
  const pred::EnergyChangePredictorSite energy_change_predictor_site_;
};

}    // namespace ansys

#endif    //LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
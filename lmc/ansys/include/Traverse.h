#ifndef LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#define LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#include "EnergyChangePredictorPairSite.h"
#include "EnergyPredictor.h"
#include "ExitTime.h"
#include "ShortRangeOrder.h"
#include "SoluteCluster.h"
#include "VacancyMigrationPredictorQuartic.h"

#include <nlohmann/json.hpp>

namespace ansys {

class Traverse {
 public:
  Traverse(unsigned long long int initial_steps,
           unsigned long long int increment_steps,
           Element solvent_element,
           std::set<Element> element_set,
           size_t smallest_cluster_criteria,
           size_t solvent_bond_criteria,
           double escape_temperature,
           const std::string &predictor_filename,
           std::string log_type,
           std::string config_type);
  virtual ~Traverse();
  void RunAnsys() const;
  void RunReformat() const;

 private:
  [[nodiscard]] std::string GetHeaderFrameString() const;
  [[nodiscard]] std::string GetHeaderClusterString() const;

  [[nodiscard]] std::string GetFrameString(const nlohmann::json &frame) const;
  [[nodiscard]] std::string GetClusterString(const nlohmann::json &frame) const;

  const unsigned long long initial_steps_;
  const unsigned long long increment_steps_;
  unsigned long long final_steps_;
  Element solvent_element_;
  const std::set<Element> element_set_;
  size_t smallest_cluster_criteria_;
  size_t solvent_bond_criteria_;
  double escape_temperature_;
  using MapVariant =
      std::variant<std::unordered_map<unsigned long long, double>, std::unordered_map<unsigned long long, std::string>>;

  std::unordered_map<std::string, MapVariant> log_map_;

  const std::string log_type_;
  const std::string config_type_;
  const pred::EnergyPredictor energy_predictor_;
  const pred::VacancyMigrationPredictorQuartic vacancy_migration_predictor_;
  const pred::EnergyChangePredictorPairSite energy_change_predictor_pair_site_;
};

}    // namespace ansys

#endif    //LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
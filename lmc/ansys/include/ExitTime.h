#ifndef LMC_LMC_ANSYS_INCLUDE_EXITTIME_H
#define LMC_LMC_ANSYS_INCLUDE_EXITTIME_H

#include "Config.h"
#include "EnergyChangePredictorPairSite.h"
#include "VacancyMigrationPredictorQuartic.h"
#include <nlohmann/json.hpp>

namespace ansys {
class ExitTime {
 public:
  ExitTime(const cfg::Config &config,
           const Element &solvent_element,
           std::set<Element> element_set,
           double temperature,
           const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor,
           const pred::EnergyChangePredictorPairSite &energy_change_predictor_site,
           const std::map<Element, double> &chemical_potential);

  void GetExitTimeInfo(nlohmann::json &frame_info,
                       std::map<std::string, cfg::Config::VectorVariant> &auxiliary_lists,
                       std::map<std::string, cfg::Config::ValueVariant> &global_list) const;

 protected:
  [[nodiscard]] double BuildMarkovChain(const std::unordered_set<size_t> &atom_id_set,
                                        const std::vector<std::vector<size_t>> &neighbor_atom_id_lists,
                                        const std::vector<std::vector<double>> &migration_barrier_lists,
                                        const std::vector<double> &base_energy_list) const;

  [[nodiscard]] std::tuple<
      std::unordered_map<std::pair<size_t, size_t>, std::pair<double, double>, boost::hash<std::pair<size_t, size_t>>>,
      std::vector<std::vector<size_t>>,
      std::vector<std::vector<double>>,
      std::vector<std::vector<double>>>
  GetJumpEnergetics(const std::unordered_set<size_t> &atom_id_set) const;
  [[nodiscard]] std::tuple<std::vector<std::vector<double>>, std::vector<double>, std::vector<double>>
  GetBarrierListAndExitTime() const;
  [[nodiscard]] std::vector<double>
  GetAverageBarriers(const std::unordered_set<size_t> &atom_id_set,
                     const std::unordered_map<std::pair<size_t, size_t>,
                                              std::pair<double, double>,
                                              boost::hash<std::pair<size_t, size_t>>> &pair_energy_map) const;
  // [[nodiscard]] double GetLocalBindingEnergy() const;
  [[nodiscard]] std::map<Element, std::vector<double>> GetBindingEnergy() const;
  [[nodiscard]] std::vector<double> GetProfileEnergy() const;

  const cfg::Config &config_;
  const Element solvent_element_;
  const std::set<Element> element_set_;
  const double beta_;
  const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor_;
  const pred::EnergyChangePredictorPairSite &energy_change_predictor_pair_site_;
  const std::map<Element, double> chemical_potential_;
};
}    // namespace ansys

#endif    //LMC_LMC_ANSYS_INCLUDE_EXITTIME_H
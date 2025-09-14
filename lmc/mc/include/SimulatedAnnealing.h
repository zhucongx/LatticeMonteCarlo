#ifndef LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_
#define LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_
#include "CanonicalMcAbstract.h"
#include "EnergyChangePredictorPair.h"

#include <mpi.h>
#include <omp.h>
#include <random>
#include <utility>

namespace mc {
class SimulatedAnnealing : public McAbstract {
 public:
  SimulatedAnnealing(const Factor_t &factors,
                     Element solvent_element,
                     const std::map<Element, size_t> &solute_atom_count,
                     unsigned long long int log_dump_steps,
                     unsigned long long int config_dump_steps,
                     unsigned long long int maximum_steps,
                     double initial_temperature,
                     const std::string &json_coefficients_filename);
  void Simulate() override;
 private:
  virtual void Dump() const;
  void UpdateTemperature(bool accepted);
  std::pair<size_t, size_t> GenerateLatticeIdJumpPair();
  bool SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair, double dE);
  // predictors / selectors
  const pred::EnergyChangePredictorPairSite energy_change_predictor_;
  mutable std::uniform_int_distribution<size_t> atom_index_selector_;

  // runtime tracking (energies / progress)
  mutable double lowest_energy_{};
  double recent_best_energy_{};
  unsigned long long int last_improvement_step_{};

  // schedule parameters
  const double initial_temperature_;

  // smooth per-step geometric: T *= exp(-baseline_c_/N)
  // choose baseline_c_ so that T_end/T0 ≈ exp(-baseline_c_)
  const double baseline_c_{3};

  // reheat ratio configuration
  const double reheat_trigger_ratio_{0.05};
  const double reheat_cooldown_ratio_{0.10};
  // acceptance window configuration
  const double window_ratio_{0.001};

  // reheat thresholds (absolute steps)
  const unsigned long long int reheat_trigger_steps_;   // steps without improvement to trigger reheat
  const unsigned long long int reheat_cooldown_steps_;  // minimal interval between reheats
  const unsigned int window_size_;                      // sliding window size for acceptance ratio

  // (min, max) targets at early and late progress
  const double acc_max_{0.50};
  const double acc_min_{0.05};
  const double cool_down_factor_{0.99};  // if acceptance too high

  // acceptance window state
  unsigned int window_trials_{0};
  unsigned int window_accepts_{0};

  // reheating configuration / state
  const double reheat_factor_{1.10};
  const size_t max_reheats_{5};
  unsigned long long int last_reheat_step_{0};
  size_t reheats_done_{0};
};
}    // namespace mc

#endif    //LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_

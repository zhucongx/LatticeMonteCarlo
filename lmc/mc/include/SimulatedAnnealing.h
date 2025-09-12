#ifndef LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_
#define LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_
#include "CanonicalMcAbstract.h"
#include "EnergyChangePredictorPair.h"

#include <mpi.h>
#include <omp.h>
#include <random>

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
  void UpdateTemperature();
  std::pair<size_t, size_t> GenerateLatticeIdJumpPair();
  bool SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair, double dE);
  // helpful properties
  const pred::EnergyChangePredictorPairSite energy_change_predictor_;
  mutable std::uniform_int_distribution<size_t> atom_index_selector_;
  mutable double lowest_energy_;
  // simulation parameters
  const double initial_temperature_;

  // internal schedule state (no external params)
  // stagewise geometric controls
  unsigned long long int stage_start_step_{0};
  const size_t num_stages_{40};
  const unsigned long long int stage_length_;
  const double stage_ratio_{0.98};
  // continuous baseline geometric cooling: T *= exp(-baseline_c_/maximum_steps_) per step
  const double baseline_c_{1.5};
  // acceptance window controls
  unsigned int window_trials_{0};
  unsigned int window_accepts_{0};
  const unsigned int window_size_{5000};
  // acceptance targets (linearly interpolated early -> late)
  const double acc_early_min_{0.35};
  const double acc_early_max_{0.55};
  const double acc_late_min_{0.05};
  const double acc_late_max_{0.20};
  const double cool_down_factor_{0.98};  // if acceptance too high
  // reheating on stagnation
  double best_energy_so_far_{0.0};
  unsigned long long int last_improvement_step_{0};
  const unsigned long long int stagnation_threshold_{};  // in steps, set relative to stage_length_
  unsigned long long int last_reheat_step_{0};
  size_t reheats_done_{0};
  const size_t max_reheats_{3};
};
}    // namespace mc

#endif    //LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_

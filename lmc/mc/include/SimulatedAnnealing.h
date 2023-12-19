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
                     double decrement_temperature,
                     const std::string &json_coefficients_filename);
  void Simulate() override;
 private:
  virtual void Dump() const;
  void UpdateTemperature();
  std::pair<size_t, size_t> GenerateLatticeIdJumpPair();
  void SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair, double dE);
  // helpful properties
  const pred::EnergyChangePredictorPairSite energy_change_predictor_;
  mutable std::uniform_int_distribution<size_t> atom_index_selector_;
  mutable double lowest_energy_;
  // simulation parameters
  const double initial_temperature_;
  const double decrement_temperature_;
};
}    // namespace mc

#endif    //LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_

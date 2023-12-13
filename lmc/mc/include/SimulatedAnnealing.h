#ifndef LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_
#define LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_
#include "CanonicalMcAbstract.h"
#include "EnergyChangePredictorPair.h"

#include <mpi.h>
#include <omp.h>
#include <random>

namespace mc {
class SimulatedAnnealing : private CanonicalMcAbstract {
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
  void Dump() const override;
  void UpdateTemperature();

  // simulation parameters
  const double initial_temperature_;
  const double decrement_temperature_;
};
}    // namespace mc

#endif    //LMC_LMC_MC_INCLUDE_SIMULATEDANNEALING_H_

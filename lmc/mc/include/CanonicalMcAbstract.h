#ifndef LMC_LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#include "EnergyChangePredictorPairSite.h"
#include "McAbstract.h"

#include <mpi.h>
#include <omp.h>
#include <random>

namespace mc {

class CanonicalMcAbstract : public McAbstract {
 public:
  CanonicalMcAbstract(cfg::Config config,
                      unsigned long long int log_dump_steps,
                      unsigned long long int config_dump_steps,
                      unsigned long long int maximum_steps,
                      unsigned long long int thermodynamic_averaging_steps,
                      unsigned long long int restart_steps,
                      double restart_energy,
                      double temperature,
                      const std::set<Element> &element_set,
                      const std::string &json_coefficients_filename);
  void Simulate() override = 0;

 protected:
  virtual void Dump() const;
  std::pair<size_t, size_t> GenerateLatticeIdJumpPair();
  void SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair, double dE);

  // helpful properties
  const pred::EnergyChangePredictorPairSite energy_change_predictor_;
  std::vector<size_t> solute_lattice_id_vector_;

  mutable std::uniform_int_distribution<size_t> atom_index_selector_;

};

}    // namespace mc

#endif    //LMC_LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_

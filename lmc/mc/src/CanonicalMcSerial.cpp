#include "CanonicalMcSerial.h"
#include <utility>
#include <chrono>
#include <omp.h>
namespace mc {
CanonicalMcSerial::CanonicalMcSerial(cfg::Config config,
                                     const unsigned long long int log_dump_steps,
                                     const unsigned long long int config_dump_steps,
                                     const unsigned long long int maximum_steps,
                                     const unsigned long long int thermodynamic_averaging_steps,
                                     const unsigned long long int restart_steps,
                                     const double restart_energy,
                                     const double temperature,

    // const double initial_temperature,
    // const double decrement_temperature,
                                     const std::set<Element> &element_set,
                                     const std::string &json_coefficients_filename)
    : CanonicalMcAbstract(std::move(config),
                          log_dump_steps,
                          config_dump_steps,
                          maximum_steps,
                          thermodynamic_averaging_steps,
                          restart_steps,
                          restart_energy,
                          temperature,
    // initial_temperature,
    // decrement_temperature,
                          element_set,
                          json_coefficients_filename) {
  if (world_size_ != 1) {
    std::cout << "Must use 1 precesses. Terminating...\n" << std::endl;
    MPI_Finalize();
    exit(0);
  }
#pragma omp parallel  default(none) shared(std::cout)
  {
#pragma omp master
    {
      std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
    }
  }
}
void CanonicalMcSerial::Simulate() {
  // while (steps_ <= maximum_steps_ * static_cast<unsigned long long int>(initial_temperature_ / decrement_temperature_ + 1)) {
  while (steps_ <= maximum_steps_) {
    auto lattice_id_jump_pair = GenerateLatticeIdJumpPair();
    auto dE = energy_change_predictor_.GetDeFromLatticeIdPair(config_, lattice_id_jump_pair);
    thermodynamic_averaging_.AddEnergy(energy_);
    Dump();
    // UpdateTemperature();
    SelectEvent(lattice_id_jump_pair, dE);
    ++steps_;
  }
}

} // mc

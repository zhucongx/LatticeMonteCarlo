#ifndef LKMC_LKMC_ANSYS_INCLUDE_SIMULATEDANNEALING_H_
#define LKMC_LKMC_ANSYS_INCLUDE_SIMULATEDANNEALING_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "StateChangePredictor.h"
namespace ansys {
constexpr double kBoltzmannConstant = 8.617333262145e-5;

class SimulatedAnnealing {
  public:
    SimulatedAnnealing(const Factor_t &factors,
                       Element solvent_element,
                       const std::map<Element, size_t> &solute_atom_count,
                       unsigned long long int log_dump_steps,
                       unsigned long long int config_dump_steps,
                       unsigned long long int maximum_number,
                       unsigned long long int early_stop_number,
                       const std::string &json_coefficients_filename);
    void Simulate();
  private:
    double MoveOneStep();
    inline void Dump(std::ofstream &ofs);
    bool EarlyStop();
    std::pair<size_t, size_t> GenerateAtomIdJumpPair();
    // simulation parameters
    cfg::Config config_;
    const std::vector<size_t> solute_atom_id_vector_;
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_number_;
    const unsigned long long int early_stop_number_;
    // simulation statistics
    unsigned long long int steps_;
    unsigned long long int count_;
    double energy_{0.0};
    double lowest_energy_{0.0};
    const double initial_temperature_{1000 * kBoltzmannConstant * 7};
    double temperature_{};

    const pred::StateChangePredictor energy_predictor_;
    mutable std::mt19937_64 generator_;
    mutable std::uniform_int_distribution<size_t> solute_atom_selector_;
    mutable std::uniform_int_distribution<size_t> neighbor_index_selector_;
    mutable std::uniform_real_distribution<double> one_distribution_{0.0, 1.0};
};
} // namespace ansys

#endif //LKMC_LKMC_ANSYS_INCLUDE_SIMULATEDANNEALING_H_

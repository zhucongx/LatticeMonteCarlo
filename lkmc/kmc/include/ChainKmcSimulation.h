#ifndef LKMC_LKMC_KMC_INCLUDE_CHAINKMCSIMULATION_H_
#define LKMC_LKMC_KMC_INCLUDE_CHAINKMCSIMULATION_H_
#include <random>
#include <mpi.h>

#include "EnergyPredictorE0DE.h"
#include "EnergyPredictorE0DEBond.h"
#include "KmcEvent.h"
namespace kmc {
//  j -> k -> i ->l
//       |
// current position
class ChainKMCSimulation {
  public:
    ChainKMCSimulation(cfg::Config config,
                       unsigned long long int log_dump_steps,
                       unsigned long long int config_dump_steps,
                       unsigned long long int maximum_number,
                       const std::set<Element> &type_set,
                       unsigned long long int steps,
                       double energy,
                       double time,
                       const std::string &json_parameters_filename);
    virtual ~ChainKMCSimulation();
    virtual void Simulate();

  protected:
    virtual bool CheckAndSolveEquilibrium(std::ofstream &ofs) { return false; };
    inline void Dump(std::ofstream &ofs);
    KMCEvent GetEventI();
    [[nodiscard]] double BuildEventListParallel();

    std::vector<size_t> GetLIndexes();
    size_t SelectEvent() const;

    // constants
    static constexpr size_t kFirstEventListSize = constants::kNumFirstNearestNeighbors;
    static constexpr size_t kSecondEventListSize = 11;

    // simulation parameters
    cfg::Config config_;
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_number_;

    // simulation statistics
    unsigned long long int steps_;
    double energy_;
    double time_;
    const size_t vacancy_index_;

    // simulation variables
    double one_step_barrier_{0.0};
    double one_step_energy_change_{0.0};
    double one_step_time_change_{0.0};

    // helpful properties
    double total_rate_k_{0.0}, total_rate_i_{0.0};
    int world_rank_{-1}, first_group_rank_{-1}, second_group_rank_{-1};
    // double rij{0}, pij{0};
    std::pair<size_t, size_t> atom_id_jump_pair_;
    size_t previous_j_;

    MPI_Group world_group_, first_group_, second_group_;
    MPI_Comm first_comm_, second_comm_;

    std::vector<KMCEvent> event_list_{};
    const pred::EnergyPredictorE0DEBond energy_predictor_;
    mutable std::mt19937_64 generator_;
};

} // namespace kmc


#endif //LKMC_LKMC_KMC_INCLUDE_CHAINKMCSIMULATION_H_

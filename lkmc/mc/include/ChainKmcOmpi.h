#ifndef LKMC_LKMC_KMC_INCLUDE_CHAINKMCOMPI_H_
#define LKMC_LKMC_KMC_INCLUDE_CHAINKMCOMPI_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "VacancyMigrationPredictorQuarticLru.h"
#include "JumpEvent.h"
namespace mc {
//  j -> k -> i -> l
//       |
// current position

class ChainKmcOmpi {
  public:
    ChainKmcOmpi(cfg::Config config,
                 unsigned long long int log_dump_steps,
                 unsigned long long int config_dump_steps,
                 unsigned long long int maximum_number,
                 double temperature,
                 const std::set<Element> &element_set,
                 unsigned long long int restart_steps,
                 double restart_energy,
                 double restart_time,
                 const std::string &json_coefficients_filename);
    virtual ~ChainKmcOmpi();
    virtual void Simulate();
  protected:
    inline void Dump(std::ofstream &ofs) const;
    void BuildFirstEventKIAndGetTotalRates();
    double UpdateIndirectProbabilityAndCalculateTime();
    size_t SelectEvent() const;
    // constants
    static constexpr size_t kEventListSize = constants::kNumFirstNearestNeighbors;

    // config
    cfg::Config config_;

    // simulation parameters
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_number_;
    double beta_;
    // simulation statistics
    unsigned long long int steps_;
    double energy_;
    double time_;
    const size_t vacancy_index_;

    // simulation variables
    double one_step_barrier_{0.0};
    double one_step_energy_change_{0.0};
    double one_step_time_change_{0.0};
    Element migrating_element_{};

    // helpful properties
    double total_rate_k_{0.0}; // k would be same for all
    double total_rate_i_{0.0}; // i would be different
    std::array<size_t, kEventListSize> l_index_list_{};
    JumpEvent event_k_i_; // first event, different for each process
    std::vector<JumpEvent> event_k_i_list_{};

    std::pair<size_t, size_t> atom_id_jump_pair_;
    size_t previous_j_;

    const pred::VacancyMigrationPredictorQuarticLru energy_predictor_;
    mutable std::mt19937_64 generator_;

  private:
    int world_rank_{-1};
    MPI_Op mpi_op_{};
    MPI_Datatype mpi_datatype_{};
};

} // namespace mc

#endif //LKMC_LKMC_KMC_INCLUDE_CHAINKMCOMPI_H_

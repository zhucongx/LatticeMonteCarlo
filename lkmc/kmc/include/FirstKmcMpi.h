#ifndef LKMC_LKMC_KMC_INCLUDE_FIRSTKMCMPI_H_
#define LKMC_LKMC_KMC_INCLUDE_FIRSTKMCMPI_H_
#include <random>
#include <mpi.h>
#include "VacancyMigrationPredictorQuarticLru.h"
#include "JumpEvent.h"
namespace kmc {
class FirstKmcMpi {
  public:
    FirstKmcMpi(cfg::Config config,
                unsigned long long int log_dump_steps,
                unsigned long long int config_dump_steps,
                unsigned long long int maximum_number,
                double temperature,
                const std::set<Element> &element_set,
                unsigned long long int restart_steps,
                double restart_energy,
                double restart_time,
                const std::string &json_coefficients_filename);
    virtual ~FirstKmcMpi();
    virtual void Simulate();
  protected:
    inline void Dump(std::ofstream &ofs) const;
    void BuildEventListParallel();
    [[nodiscard]] size_t SelectEvent() const;

    // constants
    static constexpr size_t kEventListSize = constants::kNumFirstNearestNeighbors;

    // simulation parameters
    cfg::Config config_;
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_number_;
    const double beta_;

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
    double total_rate_{0.0};
    int world_rank_{-1};
    std::pair<size_t, size_t> atom_id_jump_pair_;

    std::vector<JumpEvent> event_list_{};
    const pred::VacancyMigrationPredictorQuarticLru energy_predictor_;
    mutable std::mt19937_64 generator_;
};
} // namespace kmc

#endif //LKMC_LKMC_KMC_INCLUDE_FIRSTKMCMPI_H_

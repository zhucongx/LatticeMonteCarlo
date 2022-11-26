#ifndef LKMC_LKMC_MC_INCLUDE_CHAINKMCOPENMP_H_
#define LKMC_LKMC_MC_INCLUDE_CHAINKMCOPENMP_H_
#include <random>
#include <omp.h>
#include "VacancyMigrationPredictorQuarticLru.h"
#include "JumpEvent.h"
namespace mc {
//  j -> k -> i ->l
//       |
// current position
class KineticMcChainOmp {
  public:
    KineticMcChainOmp(const cfg::Config &config,
                      unsigned long long int log_dump_steps,
                      unsigned long long int config_dump_steps,
                      unsigned long long int maximum_steps,
                      double temperature,
                      const std::set<Element> &element_set,
                      unsigned long long int restart_steps,
                      double restart_energy,
                      double restart_time,
                      const std::string &json_coefficients_filename);
    virtual ~KineticMcChainOmp();
    virtual void Simulate();

  protected:
    inline void Dump(std::ofstream &ofs) const;
    void BuildEventKIListAndLIndexList();
    void BuildTotalRateIList();
    double CalculateTime();
    size_t SelectEvent() const;

    // constants
    static constexpr size_t kFirstEventListSize = constants::kNumFirstNearestNeighbors;
    static constexpr size_t kSecondEventListSize = kFirstEventListSize - 1;

    // config
    std::array<cfg::Config, kFirstEventListSize * kSecondEventListSize> config_list_;

    // simulation parameters
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_steps_;
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
    double total_rate_k_{0.0};
    std::array<double, kFirstEventListSize> total_rate_i_list_{};
    std::array<size_t, kFirstEventListSize * kSecondEventListSize> l_index_list_{};
    std::array<JumpEvent, kFirstEventListSize> event_k_i_list_{};

    std::pair<size_t, size_t> atom_id_jump_pair_;
    size_t previous_j_;


    const pred::VacancyMigrationPredictorQuarticLru energy_predictor_;
    mutable std::mt19937_64 generator_;
};

} // mc

#endif //LKMC_LKMC_MC_INCLUDE_CHAINKMCOPENMP_H_

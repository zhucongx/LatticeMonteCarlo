#ifndef LKMC_LKMC_KMC_INCLUDE_CHAINKMCOPENMP_H_
#define LKMC_LKMC_KMC_INCLUDE_CHAINKMCOPENMP_H_
#include <random>
#include <omp.h>

#include "EnergyPredictorStateLru.h"
#include "JumpEvent.h"
namespace kmc {
//  j -> k -> i ->l
//       |
// current position
class ChainKmcOmp {
  public:

    static constexpr double kBoltzmannConstant = 8.617333262145e-5;
    static constexpr double kPrefactor = 1e13;

    ChainKmcOmp(const cfg::Config &config,
                unsigned long long int log_dump_steps,
                unsigned long long int config_dump_steps,
                unsigned long long int maximum_number,
                double temperature,
                const std::set<Element> &type_set,
                unsigned long long int steps,
                double energy,
                double time,
                const std::string &json_parameters_filename);
    virtual ~ChainKmcOmp();
    virtual void Simulate();

  protected:
    inline void Dump(std::ofstream &ofs);
    void BuildFirstEventList();
    void BuildSecondEventList();
    double CalculateTime();
    size_t SelectEvent() const;

    // constants
    static constexpr size_t kFirstEventListSize = 12;
    static constexpr size_t kSecondEventListSize = 11;

    // simulation parameters

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

    // helpful properties
    std::array<cfg::Config, kFirstEventListSize * kSecondEventListSize> config_list_;
    std::array<double, kFirstEventListSize> total_rate_i_list_{};
    std::array<size_t, kFirstEventListSize * kSecondEventListSize> l_index_list_{};
    double total_rate_k_{0.0};
    size_t previous_j_;

    std::array<JumpEvent, kFirstEventListSize> first_event_list_{};
    const pred::EnergyPredictorStateLru energy_predictor_;
    mutable std::mt19937_64 generator_;
};

} // namespace kmc

#endif //LKMC_LKMC_KMC_INCLUDE_CHAINKMCOPENMP_H_

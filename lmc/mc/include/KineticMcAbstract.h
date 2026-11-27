#ifndef LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#include <random>
#include <string>
#include <omp.h>
#include <mpi.h>
#include "McAbstract.h"
#include "JumpEvent.h"
#include "VacancyMigrationPredictorQuarticLru.h"
#include "VacancyMigrationPredictorE0Lru.h"
#include "TimeTemperatureInterpolator.h"
namespace mc {

class KineticMcFirstAbstract : public McAbstract {
    // it is lattice id jump pair based
  public:
    KineticMcFirstAbstract(cfg::Config config,
                           unsigned long long int log_dump_steps,
                           unsigned long long int config_dump_steps,
                           unsigned long long int maximum_steps,
                           unsigned long long int thermodynamic_averaging_steps,
                           unsigned long long int restart_steps,
                           double restart_energy,
                           double restart_time,
                           double temperature,
                           const std::set<Element> &element_set,
                           const std::string &json_coefficients_filename,
                           const std::string &time_temperature_filename,
                           bool is_early_stop,
                           bool is_solute_disp,
                           unsigned long long int seed);
    ~KineticMcFirstAbstract() override;
    KineticMcFirstAbstract(const KineticMcFirstAbstract &) = delete;
    void operator=(const mc::KineticMcFirstAbstract &) = delete;
    void Simulate() override;
  protected:
    void UpdateTemperature();
    virtual void Dump() const;
    size_t SelectEvent() const;
    void Debug(double one_step_time) const;
    virtual void BuildEventList() = 0;
    virtual double CalculateTime() = 0;
    virtual void OneStepSimulation();
    void IsEscaped();
    // Counts the selected jump if its barrier is negative forward or backward, and writes a warning when the count
    // reaches 1, 10, 100 and so on.
    void CountNegativeBarrier();
    [[nodiscard]] std::string DescribeNegativeBarriers() const;
    // constants
    static constexpr size_t kEventListSize = constants::kNumFirstNearestNeighbors;

    // helpful properties
    const pred::VacancyMigrationPredictorQuarticLru vacancy_migration_predictor_lru_;
    const pred::TimeTemperatureInterpolator time_temperature_interpolator_;
    const bool is_time_temperature_interpolator_;
    size_t vacancy_lattice_id_;
    Vector_d solute_com_trajectory_{};
    const bool is_early_stop_;
    const bool is_solute_disp_;
    std::array<JumpEvent, kEventListSize> event_k_i_list_{};
    JumpEvent event_k_i_{};
    double total_rate_k_{0.0}; // k would be same for all
    const Element solvent_element_;
    const double total_solute_mass_;
    // The selected jumps of this run, kept on rank 0: all of them, those with a negative barrier forward and those
    // with a negative barrier only backward.
    unsigned long long int num_selected_jumps_{0};
    unsigned long long int num_negative_forward_barriers_{0};
    unsigned long long int num_negative_backward_barriers_{0};
    // The lowest barrier of these jumps in either direction.
    double lowest_barrier_{0.0};
    unsigned long long int next_warning_of_negative_barriers_{1};
};
//  j -> k -> i -> l
//       |
// current position
// The sums over the twelve first jumps k -> i that a second-order step is made of. j is the site the vacancy came from.
struct SecondOrderSums {
  double beta_bar_k{0.0};
  double beta_k{0.0};
  double gamma_bar_k_j{0.0};
  double gamma_k_j{0.0};
  double beta_k_j{0.0};
  double alpha_k_j{0.0};
  double ts_numerator{0.0};
  double ts_j_numerator{0.0};
  SecondOrderSums &operator+=(const SecondOrderSums &summands);
};
class KineticMcSecondAbstract : public KineticMcFirstAbstract {
  public:
    KineticMcSecondAbstract(cfg::Config config,
                            unsigned long long int log_dump_steps,
                            unsigned long long int config_dump_steps,
                            unsigned long long int maximum_steps,
                            unsigned long long int thermodynamic_averaging_steps,
                            unsigned long long int restart_steps,
                            double restart_energy,
                            double restart_time,
                            double temperature,
                            const std::set<Element> &element_set,
                            const std::string &json_coefficients_filename,
                            const std::string &time_temperature_filename,
                            bool is_early_stop,
                            bool is_solute_disp,
                            unsigned long long int seed);
    ~KineticMcSecondAbstract() override;
    KineticMcSecondAbstract(const KineticMcSecondAbstract &) = delete;
    void operator=(const mc::KineticMcSecondAbstract &) = delete;
  protected:
    void OneStepSimulation() override;
    // What the first jump k -> i of event_k_i adds to the sums of a step. total_rate_i is the sum of the rates of the
    // twelve jumps i -> l that can follow it.
    [[nodiscard]] SecondOrderSums GetSecondOrderSummands(const JumpEvent &event_k_i, double total_rate_i) const;
    // The probability that the step ends on site i, from the summands of the jump k -> i and the sums over all twelve.
    [[nodiscard]] double GetSecondOrderProbability(const JumpEvent &event_k_i,
                                                   const SecondOrderSums &summands,
                                                   const SecondOrderSums &sums) const;
    // The mean time of the step.
    [[nodiscard]] double GetSecondOrderTime(const SecondOrderSums &sums) const;
    // helpful properties
    size_t previous_j_lattice_id_;
};
} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_

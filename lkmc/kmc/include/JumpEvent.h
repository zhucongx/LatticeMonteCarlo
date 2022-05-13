#ifndef LKMC_LKMC_KMC_INCLUDE_JUMPEVENT_H_
#define LKMC_LKMC_KMC_INCLUDE_JUMPEVENT_H_
#include <cstddef>
#include <utility>
namespace kmc {
constexpr double kBoltzmannConstant = 8.617333262145e-5;
constexpr double kPrefactor = 1e13;

class JumpEvent {
  public:
    /// Constructor
    JumpEvent();
    JumpEvent(std::pair<size_t, size_t> atom_id_jump_pair,
              const std::pair<double, double> &barrier_and_diff,
              double beta);
    /// Getter
    [[nodiscard]] const std::pair<size_t, size_t> &GetAtomIdJumpPair() const;
    [[nodiscard]] double GetForwardBarrier() const;
    [[nodiscard]] double GetForwardRate() const;
    [[nodiscard]] double GetBackwardBarrier() const;
    [[nodiscard]] double GetBackwardRate() const;
    [[nodiscard]] double GetEnergyChange() const;
    [[nodiscard]] double GetProbability() const;
    [[nodiscard]] double GetCumulativeProvability() const;
    /// Setter
    void SetProbability(double probability);
    void SetCumulativeProbability(double cumulative_probability);
    void CalculateProbability(double total_rates);
  private:
    double beta_{};
    std::pair<size_t, size_t> atom_id_jump_pair_{};
    double barrier_{};
    double energy_change_{};
    double forward_rate_{};

    double probability_{};
    double cumulative_probability_{};
};
} // namespace kmc
#endif //LKMC_LKMC_KMC_INCLUDE_JUMPEVENT_H_

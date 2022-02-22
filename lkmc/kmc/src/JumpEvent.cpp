#include <cmath>
#include <utility>
#include "JumpEvent.h"
namespace kmc {
JumpEvent::JumpEvent() = default;
JumpEvent::JumpEvent(std::pair<size_t, size_t> atom_id_jump_pair,
                     const std::pair<double, double> &barrier_and_diff,
                     double coefficient)
    : coefficient_(coefficient),
      atom_id_jump_pair_(std::move(atom_id_jump_pair)),
      barrier_(barrier_and_diff.first),
      energy_change_(barrier_and_diff.second),
      forward_rate_(std::exp(-barrier_ * coefficient_)) {}
const std::pair<size_t, size_t> &JumpEvent::GetAtomIdJumpPair() const {
  return atom_id_jump_pair_;
}
double JumpEvent::GetForwardBarrier() const {
  return barrier_;
}
double JumpEvent::GetForwardRate() const {
  return forward_rate_;
}
double JumpEvent::GetBackwardBarrier() const {
  return barrier_ - energy_change_;
}
double JumpEvent::GetBackwardRate() const {
  return std::exp((energy_change_ - barrier_) * coefficient_);
}
double JumpEvent::GetEnergyChange() const {
  return energy_change_;
}

double JumpEvent::GetProbability() const {
  return probability_;
}

double JumpEvent::GetCumulativeProvability() const {
  return cumulative_probability_;
}

void JumpEvent::SetProbability(double probability) {
  probability_ = probability;
}

void JumpEvent::SetCumulativeProvability(double cumulative_provability) {
  cumulative_probability_ = cumulative_provability;
}

void JumpEvent::CalculateProbability(double total_rates) {
  probability_ = forward_rate_ / total_rates;
}
} // namespace kmc

#include <cmath>
#include <utility>
#include "JumpEvent.h"
namespace mc {
JumpEvent::JumpEvent() = default;
JumpEvent::JumpEvent(std::pair<size_t, size_t> jump_pair,
                     const std::pair<double, double> &barrier_and_diff,
                     double beta)
    : beta_(beta),
      jump_pair_(std::move(jump_pair)),
      barrier_(barrier_and_diff.first),
      energy_change_(barrier_and_diff.second),
      forward_rate_(std::exp(-barrier_ * beta_)) {}
const std::pair<size_t, size_t> &JumpEvent::GetIdJumpPair() const {
  return jump_pair_;
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
  return std::exp((energy_change_ - barrier_) * beta_);
}
double JumpEvent::GetEnergyChange() const {
  return energy_change_;
}

double JumpEvent::GetProbability() const {
  return probability_;
}

double JumpEvent::GetCumulativeProbability() const {
  return cumulative_probability_;
}

void JumpEvent::SetProbability(double probability) {
  probability_ = probability;
}

void JumpEvent::SetCumulativeProbability(double cumulative_probability) {
  cumulative_probability_ = cumulative_probability;
}

void JumpEvent::CalculateProbability(double total_rates) {
  probability_ = forward_rate_ / total_rates;
}
JumpEvent JumpEvent::GetReverseJumpEvent() const {
  return JumpEvent{{jump_pair_.second, jump_pair_.first},
                   {barrier_ - energy_change_, -energy_change_}, beta_};
}
} // mc

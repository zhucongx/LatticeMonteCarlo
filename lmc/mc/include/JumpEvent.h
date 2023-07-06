/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 4:05 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_MC_INCLUDE_JUMPEVENT_H_
#define LMC_MC_INCLUDE_JUMPEVENT_H_
#include <cstddef>
#include <utility>
namespace mc {

class JumpEvent {
 public:
  /// Constructor
  JumpEvent();
  JumpEvent(std::pair<size_t, size_t> jump_pair,
            const std::pair<double, double> &barrier_and_diff,
            double beta);
  /// Getter
  [[nodiscard]] const std::pair<size_t, size_t> &GetIdJumpPair() const;
  [[nodiscard]] double GetForwardBarrier() const;
  [[nodiscard]] double GetForwardRate() const;
  [[nodiscard]] double GetBackwardBarrier() const;
  [[nodiscard]] double GetBackwardRate() const;
  [[nodiscard]] double GetEnergyChange() const;
  [[nodiscard]] double GetProbability() const;
  [[nodiscard]] double GetCumulativeProvability() const;
  [[nodiscard]] JumpEvent GetReverseJumpEvent() const;
  /// Setter
  void SetProbability(double probability);
  void SetCumulativeProbability(double cumulative_probability);
  void CalculateProbability(double total_rates);
 private:
  double beta_{};
  std::pair<size_t, size_t> jump_pair_{};
  double barrier_{};
  double energy_change_{};
  double forward_rate_{};

  double probability_{};
  double cumulative_probability_{};
};
} // mc
#endif //LMC_MC_INCLUDE_JUMPEVENT_H_

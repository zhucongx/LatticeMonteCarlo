#ifndef LMC_LMC_PRED_INCLUDE_RATECORRECTOR_HPP_
#define LMC_LMC_PRED_INCLUDE_RATECORRECTOR_HPP_
#include <cmath>
#include <utility>
#include <vector>
#include <iostream>

#include "Constants.hpp"
namespace pred {

class RateCorrector {

  public:
    RateCorrector(double vacancy_concentration, double solute_concentration)
        : vacancy_concentration_(vacancy_concentration),
          solute_concentration_(solute_concentration) {}
    [[nodiscard]] inline double GetTimeCorrectionFactor(double temperature) const {
      return vacancy_concentration_ / GetCorrectVacancyConcentration(temperature)
          / (1 - 13 * solute_concentration_);
    }
  private:
    static inline double GetCorrectVacancyConcentration(double temperature) {
      return 1.64 * std::exp(-(0.66 / constants::kBoltzmann / temperature - 0.7));
    }
    double vacancy_concentration_;
    double solute_concentration_;
};

} // pred

#endif //LMC_LMC_PRED_INCLUDE_RATECORRECTOR_HPP_

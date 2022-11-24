#include "ThermodynamicAveraging.h"
#include <cmath>
#include "Constants.hpp"
namespace mc {
mc::ThermodynamicAveraging::ThermodynamicAveraging(size_t size)
    : size_(size) {
}
void ThermodynamicAveraging::AddEnergy(double value) {
  energy_list_.push_back(value);
  if (energy_list_.size() > size_) {
    energy_list_.pop_front();
  }
}
double ThermodynamicAveraging::GetAverage() const {
  if (energy_list_.empty()) {
    return 0;
  }
  double sum = 0.0;
  for (auto energy: energy_list_) {
    sum += energy;
  }
  return sum / static_cast<double>(energy_list_.size());
}
double ThermodynamicAveraging::GetThermodynamicAverage(double temperature) const {
  const auto average = GetAverage();
  double partition = 0.0;
  double thermodynamic_average_energy = 0.0;
  for (auto energy: energy_list_) {
    energy -= average;
    const double exp_value = std::exp(-energy / (constants::kBoltzmannConstant * temperature));
    thermodynamic_average_energy += energy * exp_value;
    partition += exp_value;
  }
  return thermodynamic_average_energy / partition + average;
}
} // mc

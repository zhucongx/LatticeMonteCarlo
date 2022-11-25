#include "ThermodynamicAveraging.h"
#include <cmath>
#include "Constants.hpp"
namespace mc {
mc::ThermodynamicAveraging::ThermodynamicAveraging(size_t size)
    : size_(size) {
}
void ThermodynamicAveraging::AddEnergy(double value) {
  if (energy_list_.size() == size_) {
    sum_ -= energy_list_.front();
    energy_list_.pop_front();
  }
  energy_list_.push_back(value);
  sum_ += value;
}
double ThermodynamicAveraging::GetAverage() const {
  if (energy_list_.empty()) {
    return 0;
  }
  return sum_ / static_cast<double>(energy_list_.size());
}
double ThermodynamicAveraging::GetThermodynamicAverage(double temperature) const {
  const auto average = GetAverage();
  double partition = 0.0;
  double thermodynamic_average_energy = 0.0;
  double beta = 1.0 / (constants::kBoltzmann * temperature);
  for (auto energy: energy_list_) {
    energy -= average;
    const double exp_value = std::exp(-energy * beta);
    thermodynamic_average_energy += energy * exp_value;
    partition += exp_value;
  }
  return thermodynamic_average_energy / partition + average;
}
} // mc

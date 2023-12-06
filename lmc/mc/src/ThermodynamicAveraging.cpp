#include "ThermodynamicAveraging.h"
#include "Constants.hpp"
#include <cmath>
namespace mc {
mc::ThermodynamicAveraging::ThermodynamicAveraging(size_t size) : size_(size) {}
void ThermodynamicAveraging::AddEnergy(const double value)
{
  if (size_ == 0) {
    return;
  }
  if (energy_list_.size() == size_) {
    sum_ -= energy_list_.front();
    energy_list_.pop_front();
  }
  energy_list_.push_back(value);
  sum_ += value;
}
double ThermodynamicAveraging::GetAverage() const
{
  if (energy_list_.empty()) {
    return 0;
  }
  return sum_ / static_cast<double>(energy_list_.size());
}
double ThermodynamicAveraging::GetThermodynamicAverage(const double beta) const
{
  if (size_ == 0) {
    return 0;
  }
  const auto average = GetAverage();
  double partition = 0.0;
  double thermodynamic_average_energy = 0.0;
  for (auto energy: energy_list_) {
    energy -= average;
    const double exp_value = std::exp(-energy * beta);
    thermodynamic_average_energy += energy * exp_value;
    partition += exp_value;
  }
  return thermodynamic_average_energy / partition + average;
}
}    // namespace mc

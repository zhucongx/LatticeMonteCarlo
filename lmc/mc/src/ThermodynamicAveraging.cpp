#include "ThermodynamicAveraging.h"
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
double ThermodynamicAveraging::GetThermodynamicAverage() const
{
  // The window holds energies of states that the Monte Carlo trajectory already visits with their Boltzmann
  // weight, so the canonical average is their plain mean. Weighting the samples by exp(-beta * E) once more would
  // count that weight twice and return, in practice, the lowest energy of the window.
  return GetAverage();
}
}    // namespace mc

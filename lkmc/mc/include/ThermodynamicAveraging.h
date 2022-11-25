#ifndef LKMC_LKMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
#define LKMC_LKMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
#include <list>
#include <deque>
#include <iterator>

namespace mc {

class ThermodynamicAveraging {
  public:
    explicit ThermodynamicAveraging(size_t size);
    void AddEnergy(double value);
    [[nodiscard]] double GetThermodynamicAverage(double temperature) const;
  private:
    [[nodiscard]] double GetAverage() const;
    std::deque<double> energy_list_;
    const size_t size_;
    double sum_{};
};

} // mc

#endif //LKMC_LKMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_

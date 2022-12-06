#ifndef LMC_LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
#define LMC_LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
#include <list>
#include <deque>
#include <iterator>

namespace mc {

class ThermodynamicAveraging {
  public:
    explicit ThermodynamicAveraging(size_t size);
    void AddEnergy(double value);
    [[nodiscard]] double GetThermodynamicAverage(double beta) const;
  private:
    [[nodiscard]] double GetAverage() const;
    std::deque<double> energy_list_;
    const size_t size_;
    double sum_{};
};

} // mc

#endif //LMC_LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_

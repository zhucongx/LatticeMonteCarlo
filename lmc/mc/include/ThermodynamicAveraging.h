/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 4:05 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
#define LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
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
  std::deque<double> energy_list_{};
  const size_t size_;
  double sum_{};
};

} // mc

#endif //LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_

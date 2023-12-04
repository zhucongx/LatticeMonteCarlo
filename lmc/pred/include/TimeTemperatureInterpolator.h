/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 4:05 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_
#define LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_
#include <utility>
#include <vector>
#include <iostream>
namespace pred {

class TimeTemperatureInterpolator {
 public:
  explicit TimeTemperatureInterpolator(const std::string &time_temperature_filename);

  explicit TimeTemperatureInterpolator(const std::vector<std::pair<double, double>> &points);

  //Computes the corresponding Y value for X using linear interpolation
  [[nodiscard]] double GetTemperature(double time) const;

 private:
  void SortPoints();
  //Our container of (x,y) data points
  //std::pair::<double, double>.first = x value
  //std::pair::<double, double>.second = y value
  std::vector<std::pair<double, double> > points_{};
};

} // pred

#endif //LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_

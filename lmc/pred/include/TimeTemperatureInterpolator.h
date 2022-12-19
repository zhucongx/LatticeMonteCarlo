#ifndef LMC_LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_
#define LMC_LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_
#include <utility>
#include <vector>
namespace pred {

class TimeTemperatureInterpolator {
  public:
    explicit TimeTemperatureInterpolator(const std::string &time_temperature_filename);

    explicit TimeTemperatureInterpolator(const std::vector<std::pair<double, double>> &points);

    //Computes the corresponding Y value for X using linear interpolation
    [[nodiscard]] double GetTemperature(double x) const;

  private:
    void SortPoints();
    //Our container of (x,y) data points
    //std::pair::<double, double>.first = x value
    //std::pair::<double, double>.second = y value
    std::vector<std::pair<double, double>> points_{};
};

} // pred

#endif //LMC_LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_

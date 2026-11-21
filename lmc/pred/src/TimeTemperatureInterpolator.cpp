#include "TimeTemperatureInterpolator.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace pred {
TimeTemperatureInterpolator::TimeTemperatureInterpolator(const std::string &time_temperature_filename)
{
  if (time_temperature_filename.empty()) { return; }
  std::ifstream ifs(time_temperature_filename, std::ifstream::in);
  if (!ifs.is_open()) {
    throw std::runtime_error("Cannot open " + time_temperature_filename);
  }
  // Every line that starts with two numbers is a (time, temperature) point. Other lines, such as a header or a
  // comment, are skipped. Reading line by line keeps the last point of a file without a final newline.
  std::string line;
  while (std::getline(ifs, line)) {
    std::istringstream line_stream(line);
    double time, temperature;
    if (!(line_stream >> time >> temperature)) { continue; }
    if (temperature <= 0) {
      throw std::runtime_error("Non-positive temperature in " + time_temperature_filename + ": " + line);
    }
    points_.emplace_back(time, temperature);
  }
  if (points_.empty()) {
    throw std::runtime_error("No (time, temperature) points found in " + time_temperature_filename);
  }
  SortPoints();
}
TimeTemperatureInterpolator::TimeTemperatureInterpolator(const std::vector<std::pair<double, double>> &points)
    : points_(points)
{
  SortPoints();
}
void TimeTemperatureInterpolator::SortPoints()
{
  //Defensive programming. Assume the caller has not sorted the table ascending order
  std::sort(points_.begin(), points_.end());
  //Ensure that no 2 adjacent x values are equal, lest we try to divide by zero when we interpolate.
  constexpr double EPSILON{1.0e-12};
  for (std::size_t i = 1; i < points_.size(); ++i) {
    if (const double deltaX{std::abs(points_[i].first - points_[i - 1].first)}; deltaX < EPSILON) {
      throw std::range_error("2 adjacent x values are equal." + std::to_string(points_[i].first));
    }
  }
}
double TimeTemperatureInterpolator::GetTemperature(const double time) const
{
  //Define a lambda that returns true if the time value of a point pair is < the caller's time value
  auto less_than = [](const std::pair<double, double> &point, double x) {
    return point.first < x;
  };
  //Find the first table entry whose value is >= caller's time value
  const auto iter = std::lower_bound(points_.cbegin(), points_.cend(), time, less_than);
  //If the caller's X value is greater than the largest X value in the table, we can't interpolate.
  if (iter == points_.cend()) { return (points_.cend() - 1)->second; }
  //If the caller's X value is less than the smallest X value in the table, we can't interpolate.
  if (iter == points_.cbegin() and time <= points_.cbegin()->first) { return points_.cbegin()->second; }
  //We can interpolate!
  const double upper_x{iter->first};
  const double upper_y{iter->second};
  const double lower_x{(iter - 1)->first};
  const double lower_y{(iter - 1)->second};

  const double deltaY{upper_y - lower_y};
  const double deltaX{upper_x - lower_x};
  return lower_y + ((time - lower_x) / deltaX) * deltaY;
}

}    // namespace pred

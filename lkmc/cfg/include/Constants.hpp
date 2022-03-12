#ifndef LKMC_LKMC_CFG_INCLUDE_CONSTANTS_HPP_
#define LKMC_LKMC_CFG_INCLUDE_CONSTANTS_HPP_

#include <cmath>
namespace constants {
constexpr double kFirstNearestNeighborsCutoff = 3.5;
constexpr double kSecondNearestNeighborsCutoff = 4.8;
constexpr double kThirdNearestNeighborsCutoff = 5.3;
// length between fourth nearest neighbors is double of length of first nearest neighbors
// constexpr double kFourthNearestNeighborsCutoff = 5.9;
// constexpr double kFifthNearestNeighborsCutoff = 6.5;
// constexpr double kSixthNearestNeighborsCutoff = 7.1;
// constexpr double kSeventhNearestNeighborsCutoff = 7.6;
constexpr double kNearNeighborsCutoff = kThirdNearestNeighborsCutoff;


constexpr size_t kNumThirdNearestSetSize = 60;

constexpr size_t kNumFirstNearestNeighbors = 12;
constexpr size_t kNumSecondNearestNeighbors = 6;
constexpr size_t kNumThirdNearestNeighbors = 24;
// constexpr size_t kNumFourthNearestNeighbors = 12;
// constexpr size_t kNumFifthNearestNeighbors = 24;
// constexpr size_t kNumSixthNearestNeighbors = 8;
// constexpr size_t kNumSeventhNearestNeighbors = 48;
} // namespace constants
#endif //LKMC_LKMC_CFG_INCLUDE_CONSTANTS_HPP_

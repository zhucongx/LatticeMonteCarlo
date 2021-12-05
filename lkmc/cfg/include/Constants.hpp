#ifndef LKMC_LKMC_CFG_INCLUDE_CONSTANTS_HPP_
#define LKMC_LKMC_CFG_INCLUDE_CONSTANTS_HPP_

#include <cmath>
namespace constants {
constexpr double kFirstNearestNeighborsCutoff = 3.5;
constexpr double kSecondNearestNeighborsCutoff = 4.8;
constexpr double kThirdNearestNeighborsCutoff = 5.3;
// length between fourth nearest neighbors is double of length of first nearest neighbors
constexpr double kFourthNearestNeighborsCutoff = 5.9;
constexpr double kFifthNearestNeighborsCutoff = 6.5;
constexpr double kSixthNearestNeighborsCutoff = 7.1;
constexpr double kSeventhNearestNeighborsCutoff = 7.6;
// constexpr double kNearNeighborsCutoff = kThirdNearestNeighborsCutoff;

// For FCC the first nearest neighbor is 12
constexpr size_t kNumFirstNearestNeighbors = 12;
// For FCC the second nearest neighbor is 6
constexpr size_t kNumSecondNearestNeighbors = 6;
// For FCC the third nearest neighbor is 6
constexpr size_t kNumThirdNearestNeighbors = 24;
// For FCC the fourth nearest neighbor is 12
constexpr size_t kNumFourthNearestNeighbors = 12;
// For FCC the fifth nearest neighbor is 24
constexpr size_t kNumFifthNearestNeighbors = 24;
// For FCC the sixth nearest neighbor is 12
constexpr size_t kNumSixthNearestNeighbors = 8;
// For FCC the seventh nearest neighbor is 48
constexpr size_t kNumSeventhNearestNeighbors = 48;
} // namespace constants
#endif //LKMC_LKMC_CFG_INCLUDE_CONSTANTS_HPP_

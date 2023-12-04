#ifndef LMC_LMC_CFG_INCLUDE_CONSTANTS_HPP_
#define LMC_LMC_CFG_INCLUDE_CONSTANTS_HPP_

#include <cmath>
namespace constants {
constexpr double kLatticeConstant = 4.046;

constexpr double kFirstNearestNeighborsCutoff = 3.5;
constexpr double kSecondNearestNeighborsCutoff = 4.8;
constexpr double kThirdNearestNeighborsCutoff = 5.3;
// length between fourth nearest neighbors is double of length of first nearest neighbors
constexpr double kFourthNearestNeighborsCutoff = 5.9;
constexpr double kFifthNearestNeighborsCutoff = 6.5;
constexpr double kSixthNearestNeighborsCutoff = 7.1;
constexpr double kSeventhNearestNeighborsCutoff = 7.6;
constexpr double kNearNeighborsCutoff = kThirdNearestNeighborsCutoff;

constexpr size_t kNumThirdNearestSetSizeOfPair = 60;
constexpr size_t kNumThirdNearestSetSizeOfSite = 1 + 12 + 6 + 24;

constexpr size_t kNumFirstNearestNeighbors = 12;
constexpr size_t kNumSecondNearestNeighbors = 6;
constexpr size_t kNumThirdNearestNeighbors = 24;
constexpr size_t kNumFourthNearestNeighbors = 12;
constexpr size_t kNumFifthNearestNeighbors = 24;
constexpr size_t kNumSixthNearestNeighbors = 8;
constexpr size_t kNumSeventhNearestNeighbors = 48;

constexpr double kBoltzmann = 8.617333262145e-5;
constexpr double kPrefactor = 1e13;
} // constants
#endif //LMC_LMC_CFG_INCLUDE_CONSTANTS_HPP_

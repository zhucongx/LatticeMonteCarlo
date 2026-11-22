#include "VacancyMigrationPredictorQuartic.h"
#include "cfg/test/TemporaryPath.hpp"
#include "pred/test/SyntheticModel.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <utility>

namespace pred {
namespace {
// A jump of 2.8 Angstrom between minima whose curvatures add up to 3 eV per square Angstrom.
constexpr double kD = 2.8;
constexpr double kKs = 3.0;
// The barrier is positive in both directions as long as |dE| is below this value, 1.96 eV.
constexpr double kThreshold = kKs * kD * kD / 12;

// (27 / 16384) Ks D^2 (u + 4/3)^3 (4 - u) with u = 16 dE / (Ks D^2)
double ClosedForm(const double dE) {
  const double u = 16 * dE / (kKs * kD * kD);
  return 27.0 / 16384 * kKs * kD * kD * std::pow(u + 4.0 / 3, 3) * (4 - u);
}
}    // namespace

TEST(QuarticBarrier, SymmetricJumpHasTheBarrierOfTheDoubleWell) {
  // Without an energy change the form is a (x^2 - D^2 / 4)^2 with a = Ks / (4 D^2), whose barrier is Ks D^2 / 64.
  EXPECT_NEAR(GetQuarticBarrier(0.0, kD, kKs), kKs * kD * kD / 64, 1e-12);
}

TEST(QuarticBarrier, FollowsTheClosedFormOnBothSidesOfTheThreshold) {
  for (const double fraction: {-3.5, -1.5, -1.0, -0.999, -0.5, 0.0, 0.3, 0.999, 1.0, 1.5, 2.9, 3.5}) {
    const double dE = fraction * kThreshold;
    EXPECT_NEAR(GetQuarticBarrier(dE, kD, kKs), ClosedForm(dE), 1e-12) << "dE = " << fraction << " thresholds";
  }
}

// The threshold scales with Ks D^2, so it is small for a soft or short jump.
TEST(QuarticBarrier, IsNegativeBeyondTheThreshold) {
  for (const auto &[d, ks]: {std::pair<double, double>{kD, kKs}, {1.7, 0.76}, {2.3, 1.34}, {3.6, 0.42}, {5.5, 4.0}}) {
    const double threshold = ks * d * d / 12;
    for (const double fraction: {0.0, 0.5, 0.99, 0.999}) {
      EXPECT_GT(GetQuarticBarrier(-fraction * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
      EXPECT_GT(GetQuarticBarrier(fraction * threshold, d, ks), fraction * threshold) << "D " << d << ", Ks " << ks;
    }
    for (const double fraction: {1.001, 1.01, 1.5, 2.99}) {
      EXPECT_LT(GetQuarticBarrier(-fraction * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
      EXPECT_GT(GetQuarticBarrier(fraction * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
    }
    // From three thresholds on it is negative in both directions.
    for (const double fraction: {3.01, 5.0}) {
      EXPECT_LT(GetQuarticBarrier(-fraction * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
      EXPECT_LT(GetQuarticBarrier(fraction * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
    }
  }
}

// The rates of a jump and of the jump back are in detailed balance if their barriers differ by the energy change.
TEST(QuarticBarrier, ForwardAndBackwardBarriersDifferByTheEnergyChange) {
  for (const double fraction: {0.0, 0.1, 0.5, 0.999999, 1.000001, 1.01, 1.5, 2.9, 3.5, 5.0}) {
    const double dE = fraction * kThreshold;
    EXPECT_NEAR(GetQuarticBarrier(dE, kD, kKs) - GetQuarticBarrier(-dE, kD, kKs), dE, 1e-12)
        << "dE = " << fraction << " thresholds";
  }
}

namespace {
// Gives a test access to the three quantities the barrier is computed from.
class ExposedQuarticPredictor : public VacancyMigrationPredictorQuartic {
 public:
  using VacancyMigrationPredictorQuartic::VacancyMigrationPredictorQuartic;
  using VacancyMigrationPredictorQuartic::GetD;
  using VacancyMigrationPredictorQuartic::GetDe;
  using VacancyMigrationPredictorQuartic::GetKs;
};
}    // namespace

// The predictor on a model whose energy changes are so large that four of ten jumps have a negative barrier in one
// direction. It returns the value of the form whatever its sign, and jumping back passes the same point of the path.
TEST(QuarticBarrier, PredictorReturnsTheValueOfTheForm) {
  const test::TemporaryDirectory temporary_directory;
  auto config = MakeTestAlloy();
  const auto model_filename = WriteModel(temporary_directory.GetPath() / "model_with_large_energy_changes.json",
                                         MakeSyntheticModelWithScaledEnergies(config, TestElementSet(), 600));
  const ExposedQuarticPredictor predictor(model_filename, config, TestElementSet());
  std::mt19937_64 generator(7);
  std::uniform_int_distribution<size_t> pick(0, constants::kNumFirstNearestNeighbors - 1);
  size_t num_positive = 0;
  size_t num_negative_forward = 0;
  size_t num_negative_backward = 0;
  for (size_t step = 0; step < 40; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    const auto &neighbors = config.GetFirstNeighborsAdjacencyList()[vacancy];
    for (const size_t neighbor: neighbors) {
      const std::pair<size_t, size_t> jump{vacancy, neighbor};
      const auto [barrier, energy_change] = predictor.GetBarrierAndDiffFromLatticeIdPair(config, jump);
      ASSERT_EQ(energy_change, predictor.GetDe(config, jump)) << "step " << step;
      ASSERT_EQ(barrier, GetQuarticBarrier(energy_change, predictor.GetD(config, jump), predictor.GetKs(config, jump)))
          << "step " << step;
      ++(barrier < 0.0 ? num_negative_forward : barrier < energy_change ? num_negative_backward : num_positive);
    }
    const size_t target = neighbors[pick(generator)];
    const auto [barrier, energy_change] = predictor.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, target});
    config.LatticeJump({vacancy, target});
    const auto [backward_barrier, backward_change] =
        predictor.GetBarrierAndDiffFromLatticeIdPair(config, {target, vacancy});
    ASSERT_NEAR(backward_change, -energy_change, 1e-12) << "step " << step;
    ASSERT_NEAR(backward_barrier, barrier - energy_change, 1e-9 * std::max(1.0, std::abs(barrier))) << "step " << step;
  }
  // The model was made so that all three cases occur.
  EXPECT_GE(num_positive, 50u);
  EXPECT_GE(num_negative_forward, 50u);
  EXPECT_GE(num_negative_backward, 50u);
}
}    // namespace pred

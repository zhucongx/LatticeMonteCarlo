#include "VacancyMigrationPredictorQuartic.h"
#include "VacancyMigrationPredictorQuarticLru.h"
#include "cfg/test/TemporaryPath.hpp"
#include "pred/test/SyntheticModel.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <random>

namespace pred {
namespace {
// A jump of 2.8 Angstrom between minima whose curvatures add up to 3 eV per square Angstrom.
constexpr double kD = 2.8;
constexpr double kKs = 3.0;
// The quartic form has a barrier in both directions as long as |dE| is below this value, 1.96 eV.
constexpr double kThreshold = kKs * kD * kD / 12;
// Relative distances from the threshold, from close to it to far beyond it. Up to the distance 1 they are also used
// inside the threshold; the last three lie beyond |dE| = Ks D^2 / 4, where the form gives a negative barrier in both
// directions.
//
// For a jump that lowers the energy the threshold is where the barrier changes its sign, which is resolved almost to
// the last bit. For a jump that raises it the threshold is where the barrier equals dE, a comparison of two numbers
// of the order of eV whose difference grows with the third power of the distance: it is resolved only to about 1e-5
// of the threshold.
constexpr double kMarginsOfJumpDownhill[] = {1e-6, 1e-4, 1e-2, 0.5, 1.0, 2.5, 4.0, 9.0};
constexpr double kMarginsOfJumpUphill[] = {1e-4, 1e-2, 0.5, 1.0, 2.5, 4.0, 9.0};

// (27 / 16384) Ks D^2 (u + 4/3)^3 (4 - u) with u = 16 dE / (Ks D^2)
double ClosedForm(const double dE) {
  const double u = 16 * dE / (kKs * kD * kD);
  return 27.0 / 16384 * kKs * kD * kD * std::pow(u + 4.0 / 3, 3) * (4 - u);
}

// The barrier of the form together with the answer whether the jump has a barrier forward and backward.
bool HasBarrier(const double dE, const double D, const double Ks) {
  return HasBarrierForwardAndBackward(GetQuarticModelBarrier(dE, D, Ks), dE);
}
}    // namespace

TEST(QuarticBarrier, SymmetricJumpHasTheBarrierOfTheDoubleWell) {
  // Without an energy change the form is a (x^2 - D^2 / 4)^2 with a = Ks / (4 D^2), whose barrier is Ks D^2 / 64.
  EXPECT_NEAR(GetQuarticModelBarrier(0.0, kD, kKs), kKs * kD * kD / 64, 1e-12);
  EXPECT_TRUE(HasBarrier(0.0, kD, kKs));
}

TEST(QuarticBarrier, ModelBarrierFollowsTheClosedFormOnBothSidesOfTheThreshold) {
  for (const double fraction: {-3.5, -1.5, -1.0, -0.999, -0.5, 0.0, 0.3, 0.999, 1.0, 1.5, 2.9, 3.5}) {
    const double dE = fraction * kThreshold;
    EXPECT_NEAR(GetQuarticModelBarrier(dE, kD, kKs), ClosedForm(dE), 1e-12) << "dE = " << fraction << " thresholds";
  }
}

// Inside the threshold the form has a barrier forward and backward.
TEST(QuarticBarrier, InsideTheThresholdAJumpHasABarrierInBothDirections) {
  for (const double margin: kMarginsOfJumpDownhill) {
    if (margin > 1) {
      continue;
    }
    const double dE = -kThreshold * (1 - margin);
    EXPECT_GT(GetQuarticModelBarrier(dE, kD, kKs), 0.0) << "margin " << margin;
    EXPECT_TRUE(HasBarrier(dE, kD, kKs)) << "margin " << margin;
  }
  for (const double margin: kMarginsOfJumpUphill) {
    if (margin > 1) {
      continue;
    }
    const double dE = kThreshold * (1 - margin);
    EXPECT_GT(GetQuarticModelBarrier(dE, kD, kKs), dE) << "margin " << margin;
    EXPECT_TRUE(HasBarrier(dE, kD, kKs)) << "margin " << margin;
  }
}

// Beyond the threshold on the side of negative dE the start of the jump is not a minimum of the form. The barrier is
// negative there, which as a rate would let the jump happen faster than any attempt.
TEST(QuarticBarrier, BeyondTheThresholdAJumpThatLowersTheEnergyHasNoBarrier) {
  for (const double margin: kMarginsOfJumpDownhill) {
    const double dE = -kThreshold * (1 + margin);
    EXPECT_LT(GetQuarticModelBarrier(dE, kD, kKs), 0.0) << "margin " << margin;
    EXPECT_FALSE(HasBarrier(dE, kD, kKs)) << "margin " << margin;
  }
}

// Beyond the threshold on the side of positive dE the end of the jump is not a minimum. The barrier is below dE
// there, so that the jump back has a negative barrier; it falls with growing dE and is negative itself from
// dE = Ks D^2 / 4 on.
TEST(QuarticBarrier, BeyondTheThresholdAJumpThatRaisesTheEnergyHasNoBarrierBack) {
  for (const double margin: kMarginsOfJumpUphill) {
    const double dE = kThreshold * (1 + margin);
    EXPECT_LT(GetQuarticModelBarrier(dE, kD, kKs), dE) << "margin " << margin;
    EXPECT_FALSE(HasBarrier(dE, kD, kKs)) << "margin " << margin;
  }
  EXPECT_LT(GetQuarticModelBarrier(kThreshold * 3.5, kD, kKs), 0.0);
}

TEST(QuarticBarrier, AtTheThresholdTheBarrierVanishesInOneDirection) {
  EXPECT_NEAR(GetQuarticModelBarrier(-kThreshold, kD, kKs), 0.0, 1e-12);
  EXPECT_NEAR(GetQuarticModelBarrier(kThreshold, kD, kKs), kThreshold, 1e-12);
}

// The rates of a jump and of the jump back are in detailed balance if their barriers differ by the energy change.
// The form has this property on both sides of the threshold.
TEST(QuarticBarrier, ForwardAndBackwardBarriersDifferByTheEnergyChange) {
  for (const double fraction: {0.0, 0.1, 0.5, 0.999999, 1.000001, 1.01, 1.5, 2.9, 3.5, 5.0}) {
    const double dE = fraction * kThreshold;
    const double forward = GetQuarticModelBarrier(dE, kD, kKs);
    const double backward = GetQuarticModelBarrier(-dE, kD, kKs);
    EXPECT_NEAR(forward - backward, dE, 1e-10 * std::max(1.0, std::abs(forward))) << "dE = " << fraction
                                                                                   << " thresholds";
    // A jump and the jump back either both have a barrier in both directions or both have not. The two energy changes
    // next to the threshold are left out: there the comparison of the forward barrier with dE is not resolved.
    if (std::abs(fraction - 1) > 1e-3) {
      EXPECT_EQ(HasBarrierForwardAndBackward(forward, dE), HasBarrierForwardAndBackward(backward, -dE))
          << "dE = " << fraction << " thresholds";
    }
  }
}

// The threshold scales with Ks D^2, so a soft or short jump loses its barrier at a small energy change.
TEST(QuarticBarrier, ThresholdIsKsDSquaredOverTwelveForOtherJumpsToo) {
  for (const auto &[d, ks]: {std::pair<double, double>{1.7, 0.76}, {2.3, 1.34}, {3.6, 0.42}, {5.5, 4.0}}) {
    const double threshold = ks * d * d / 12;
    EXPECT_TRUE(HasBarrier(-0.999 * threshold, d, ks)) << "D " << d << ", Ks " << ks;
    EXPECT_FALSE(HasBarrier(-1.001 * threshold, d, ks)) << "D " << d << ", Ks " << ks;
    EXPECT_LT(GetQuarticModelBarrier(-1.001 * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
    EXPECT_TRUE(HasBarrier(0.999 * threshold, d, ks)) << "D " << d << ", Ks " << ks;
    EXPECT_FALSE(HasBarrier(1.001 * threshold, d, ks)) << "D " << d << ", Ks " << ks;
    EXPECT_LT(GetQuarticModelBarrier(1.001 * threshold, d, ks), 1.001 * threshold) << "D " << d << ", Ks " << ks;
  }
}

// A barrier that is not a number is not taken for a barrier.
TEST(QuarticBarrier, BarrierThatIsNotANumberCountsAsNoBarrier) {
  EXPECT_TRUE(std::isnan(GetQuarticModelBarrier(-0.5, 0.0, kKs)));
  EXPECT_TRUE(std::isnan(GetQuarticModelBarrier(std::nan(""), kD, kKs)));
  EXPECT_FALSE(HasBarrierForwardAndBackward(std::nan(""), -0.5));
  EXPECT_FALSE(HasBarrierForwardAndBackward(0.3, std::nan("")));
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

// The predictors on a model whose energy changes are so large that four of ten jumps have no barrier in the quartic
// form.
class JumpsWithoutBarrier : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    temporary_directory_ = std::make_unique<test::TemporaryDirectory>();
    config_ = MakeTestAlloy();
    model_filename_ = WriteModel(temporary_directory_->GetPath() / "model_with_large_energy_changes.json",
                                 MakeSyntheticModelWithScaledEnergies(config_, TestElementSet(), 600));
  }

  static void TearDownTestSuite() {
    temporary_directory_.reset();
  }

  static inline std::unique_ptr<test::TemporaryDirectory> temporary_directory_{};
  static inline cfg::Config config_{};
  static inline std::string model_filename_{};
};

// The predictor returns the barrier of the form as it is, also beyond the threshold, where the caller can tell from
// the barrier and the energy change that the jump has none.
TEST_F(JumpsWithoutBarrier, PredictorReturnsTheBarrierOfTheFormOnBothSidesOfTheThreshold) {
  const ExposedQuarticPredictor predictor(model_filename_, config_, TestElementSet());
  auto config = config_;
  std::mt19937_64 generator(7);
  std::uniform_int_distribution<size_t> pick(0, constants::kNumFirstNearestNeighbors - 1);
  size_t num_regular = 0;
  size_t num_without_forward_barrier = 0;
  size_t num_without_backward_barrier = 0;
  for (size_t step = 0; step < 40; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    for (const size_t neighbor: config.GetFirstNeighborsAdjacencyList()[vacancy]) {
      const std::pair<size_t, size_t> jump{vacancy, neighbor};
      const double dE = predictor.GetDe(config, jump);
      const double D = predictor.GetD(config, jump);
      const double Ks = predictor.GetKs(config, jump);
      const double threshold = Ks * D * D / 12;
      const auto [barrier, energy_change] = predictor.GetBarrierAndDiffFromLatticeIdPair(config, jump);
      ASSERT_EQ(energy_change, dE);
      ASSERT_EQ(barrier, GetQuarticModelBarrier(dE, D, Ks));
      // Leave out a jump that sits on the threshold to within the resolution of the comparison with dE.
      if (std::abs(std::abs(dE) - threshold) < 1e-4 * threshold) {
        continue;
      }
      if (dE < -threshold) {
        ++num_without_forward_barrier;
        ASSERT_LT(barrier, 0.0) << "step " << step;
        ASSERT_FALSE(HasBarrierForwardAndBackward(barrier, dE)) << "step " << step;
      } else if (dE > threshold) {
        ++num_without_backward_barrier;
        ASSERT_LT(barrier, dE) << "step " << step;
        ASSERT_FALSE(HasBarrierForwardAndBackward(barrier, dE)) << "step " << step;
      } else {
        ++num_regular;
        ASSERT_GT(barrier, std::max(dE, 0.0)) << "step " << step;
        ASSERT_TRUE(HasBarrierForwardAndBackward(barrier, dE)) << "step " << step;
      }
    }
    config.LatticeJump({vacancy, config.GetFirstNeighborsAdjacencyList()[vacancy][pick(generator)]});
  }
  // The model was made so that all three cases occur.
  EXPECT_GE(num_regular, 50u);
  EXPECT_GE(num_without_forward_barrier, 50u);
  EXPECT_GE(num_without_backward_barrier, 50u);
}

// Jumping back passes the same stationary point of the form: the two barriers differ by the energy change whether
// the form has a barrier or not, and the jump back has none exactly if the jump has none.
TEST_F(JumpsWithoutBarrier, ReverseJumpIsConsistentWithAndWithoutABarrier) {
  const VacancyMigrationPredictorQuartic predictor(model_filename_, config_, TestElementSet());
  auto config = config_;
  std::mt19937_64 generator(11);
  std::uniform_int_distribution<size_t> pick(0, constants::kNumFirstNearestNeighbors - 1);
  size_t num_without_barrier = 0;
  for (size_t step = 0; step < 200; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    const size_t target = config.GetFirstNeighborsAdjacencyList()[vacancy][pick(generator)];
    const auto [forward_barrier, forward_change] =
        predictor.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, target});
    config.LatticeJump({vacancy, target});
    const auto [backward_barrier, backward_change] =
        predictor.GetBarrierAndDiffFromLatticeIdPair(config, {target, vacancy});
    const double tolerance = 1e-9 * std::max(1.0, std::abs(forward_barrier));
    ASSERT_NEAR(backward_change, -forward_change, 1e-12) << "step " << step;
    ASSERT_NEAR(backward_barrier, forward_barrier - forward_change, tolerance) << "step " << step;
    // Leave out a jump whose barrier vanishes in one direction to within the tolerance.
    if (std::min(std::abs(forward_barrier), std::abs(backward_barrier)) < tolerance) {
      continue;
    }
    ASSERT_EQ(HasBarrierForwardAndBackward(forward_barrier, forward_change),
              HasBarrierForwardAndBackward(backward_barrier, backward_change))
        << "step " << step;
    if (!HasBarrierForwardAndBackward(forward_barrier, forward_change)) {
      ++num_without_barrier;
    }
  }
  EXPECT_GE(num_without_barrier, 20u);
}

TEST_F(JumpsWithoutBarrier, CachedPredictorAgreesWithThePlainPredictor) {
  const VacancyMigrationPredictorQuartic plain(model_filename_, config_, TestElementSet());
  const VacancyMigrationPredictorQuarticLru cached(model_filename_, config_, TestElementSet(), 64);
  auto config = config_;
  size_t num_without_barrier = 0;
  for (size_t step = 0; step < 20; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    for (const size_t neighbor: config.GetFirstNeighborsAdjacencyList()[vacancy]) {
      const auto expected = plain.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, neighbor});
      if (!HasBarrierForwardAndBackward(expected.first, expected.second)) {
        ++num_without_barrier;
      }
      // The first call fills the cache, the second one is answered from it.
      for (int repetition = 0; repetition < 2; ++repetition) {
        const auto actual = cached.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, neighbor});
        ASSERT_EQ(actual.first, expected.first) << "step " << step;
        ASSERT_EQ(actual.second, expected.second) << "step " << step;
      }
    }
    config.LatticeJump({vacancy, config.GetFirstNeighborsAdjacencyList()[vacancy][step % 12]});
  }
  EXPECT_GE(num_without_barrier, 50u);
}
}    // namespace pred

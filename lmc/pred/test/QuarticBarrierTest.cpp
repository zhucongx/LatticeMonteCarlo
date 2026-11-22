#include "VacancyMigrationPredictorQuartic.h"
#include "VacancyMigrationPredictorQuarticLru.h"
#include "cfg/test/TemporaryPath.hpp"
#include "pred/test/SyntheticModel.hpp"

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
// For a jump that lowers the energy the threshold is where the model barrier changes its sign, which is resolved
// almost to the last bit. For a jump that raises it the threshold is where the model barrier equals dE, a comparison
// of two numbers of the order of eV whose difference grows with the third power of the distance: it is resolved only
// to about 1e-5 of the threshold. Closer than that the two barriers agree to 1e-15 eV, so it does not matter which
// of them is returned.
constexpr double kMarginsOfJumpDownhill[] = {1e-6, 1e-4, 1e-2, 0.5, 1.0, 2.5, 4.0, 9.0};
constexpr double kMarginsOfJumpUphill[] = {1e-4, 1e-2, 0.5, 1.0, 2.5, 4.0, 9.0};

// (27 / 16384) Ks D^2 (u + 4/3)^3 (4 - u) with u = 16 dE / (Ks D^2)
double ClosedForm(const double dE) {
  const double u = 16 * dE / (kKs * kD * kD);
  return 27.0 / 16384 * kKs * kD * kD * std::pow(u + 4.0 / 3, 3) * (4 - u);
}
}    // namespace

TEST(QuarticBarrier, SymmetricJumpHasTheBarrierOfTheDoubleWell) {
  // Without an energy change the form is a (x^2 - D^2 / 4)^2 with a = Ks / (4 D^2), whose barrier is Ks D^2 / 64.
  EXPECT_NEAR(GetQuarticModelBarrier(0.0, kD, kKs), kKs * kD * kD / 64, 1e-12);
  EXPECT_EQ(GetQuarticBarrier(0.0, kD, kKs), GetQuarticModelBarrier(0.0, kD, kKs));
}

TEST(QuarticBarrier, ModelBarrierFollowsTheClosedFormOnBothSidesOfTheThreshold) {
  for (const double fraction: {-3.5, -1.5, -1.0, -0.999, -0.5, 0.0, 0.3, 0.999, 1.0, 1.5, 2.9, 3.5}) {
    const double dE = fraction * kThreshold;
    EXPECT_NEAR(GetQuarticModelBarrier(dE, kD, kKs), ClosedForm(dE), 1e-12) << "dE = " << fraction << " thresholds";
  }
}

// Inside the threshold the model has a barrier forward and backward, and the barrier is that of the model to the last
// bit: a simulation that never leaves the range of the form is not changed by the limits.
TEST(QuarticBarrier, InsideTheThresholdTheBarrierIsThatOfTheModel) {
  for (const double margin: kMarginsOfJumpDownhill) {
    if (margin > 1) {
      continue;
    }
    const double dE = -kThreshold * (1 - margin);
    EXPECT_GT(GetQuarticModelBarrier(dE, kD, kKs), 0.0) << "margin " << margin;
    EXPECT_EQ(GetQuarticBarrier(dE, kD, kKs), GetQuarticModelBarrier(dE, kD, kKs)) << "margin " << margin;
  }
  for (const double margin: kMarginsOfJumpUphill) {
    if (margin > 1) {
      continue;
    }
    const double dE = kThreshold * (1 - margin);
    EXPECT_GT(GetQuarticModelBarrier(dE, kD, kKs), dE) << "margin " << margin;
    EXPECT_EQ(GetQuarticBarrier(dE, kD, kKs), GetQuarticModelBarrier(dE, kD, kKs)) << "margin " << margin;
  }
}

// Beyond the threshold on the side of negative dE the start of the jump is not a minimum of the form. The model
// barrier is negative there, which as a rate would let the jump happen faster than any attempt.
TEST(QuarticBarrier, BeyondTheThresholdAJumpThatLowersTheEnergyHasNoBarrier) {
  for (const double margin: kMarginsOfJumpDownhill) {
    const double dE = -kThreshold * (1 + margin);
    EXPECT_LT(GetQuarticModelBarrier(dE, kD, kKs), 0.0) << "margin " << margin;
    EXPECT_EQ(GetQuarticBarrier(dE, kD, kKs), 0.0) << "margin " << margin;
    // not the negative zero, which would be printed as "-0" in the log
    EXPECT_FALSE(std::signbit(GetQuarticBarrier(dE, kD, kKs))) << "margin " << margin;
  }
}

// Beyond the threshold on the side of positive dE the end of the jump is not a minimum. The model barrier is below
// dE there, falls with growing dE and is negative from dE = Ks D^2 / 4 on: the jump back would have a negative
// barrier. Such a jump costs exactly its energy change.
TEST(QuarticBarrier, BeyondTheThresholdAJumpThatRaisesTheEnergyCostsItsEnergyChange) {
  for (const double margin: kMarginsOfJumpUphill) {
    const double dE = kThreshold * (1 + margin);
    EXPECT_LT(GetQuarticModelBarrier(dE, kD, kKs), dE) << "margin " << margin;
    EXPECT_EQ(GetQuarticBarrier(dE, kD, kKs), dE) << "margin " << margin;
  }
  EXPECT_LT(GetQuarticModelBarrier(kThreshold * 3.5, kD, kKs), 0.0);
}

TEST(QuarticBarrier, AtTheThresholdTheBarrierOfTheModelVanishes) {
  EXPECT_NEAR(GetQuarticModelBarrier(-kThreshold, kD, kKs), 0.0, 1e-12);
  EXPECT_NEAR(GetQuarticModelBarrier(kThreshold, kD, kKs), kThreshold, 1e-12);
  EXPECT_GE(GetQuarticBarrier(-kThreshold, kD, kKs), 0.0);
  EXPECT_NEAR(GetQuarticBarrier(-kThreshold, kD, kKs), 0.0, 1e-12);
  EXPECT_GE(GetQuarticBarrier(kThreshold, kD, kKs), kThreshold);
  EXPECT_NEAR(GetQuarticBarrier(kThreshold, kD, kKs), kThreshold, 1e-12);
}

// The rates of a jump and of the jump back are in detailed balance if their barriers differ by the energy change.
// The model has this property, and the limits keep it where the model has no barrier.
TEST(QuarticBarrier, ForwardAndBackwardBarriersDifferByTheEnergyChange) {
  for (const double fraction: {0.0, 0.1, 0.5, 0.999999, 1.000001, 1.01, 1.5, 2.9, 3.5, 5.0}) {
    const double dE = fraction * kThreshold;
    const double forward = GetQuarticBarrier(dE, kD, kKs);
    const double backward = GetQuarticBarrier(-dE, kD, kKs);
    EXPECT_NEAR(forward - backward, dE, 1e-12) << "dE = " << fraction << " thresholds";
    EXPECT_GE(forward, 0.0);
    EXPECT_GE(backward, 0.0);
  }
}

// The threshold scales with Ks D^2, so a soft or short jump loses its barrier at a small energy change.
TEST(QuarticBarrier, ThresholdIsKsDSquaredOverTwelveForOtherJumpsToo) {
  for (const auto &[d, ks]: {std::pair<double, double>{1.7, 0.76}, {2.3, 1.34}, {3.6, 0.42}, {5.5, 4.0}}) {
    const double threshold = ks * d * d / 12;
    EXPECT_GT(GetQuarticModelBarrier(-0.999 * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
    EXPECT_LT(GetQuarticModelBarrier(-1.001 * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
    EXPECT_EQ(GetQuarticBarrier(-1.001 * threshold, d, ks), 0.0) << "D " << d << ", Ks " << ks;
    EXPECT_GT(GetQuarticModelBarrier(0.999 * threshold, d, ks), 0.999 * threshold) << "D " << d << ", Ks " << ks;
    EXPECT_LT(GetQuarticModelBarrier(1.001 * threshold, d, ks), 1.001 * threshold) << "D " << d << ", Ks " << ks;
    EXPECT_EQ(GetQuarticBarrier(1.001 * threshold, d, ks), 1.001 * threshold) << "D " << d << ", Ks " << ks;
  }
}

// A barrier that is not a number has to stay visible instead of being replaced by one of the limits.
TEST(QuarticBarrier, BarrierThatIsNotANumberIsNotHidden) {
  EXPECT_TRUE(std::isnan(GetQuarticBarrier(-0.5, 0.0, kKs)));
  EXPECT_TRUE(std::isnan(GetQuarticBarrier(std::nan(""), kD, kKs)));
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

TEST_F(JumpsWithoutBarrier, PredictorLimitsTheBarrierExactlyWhereTheModelHasNone) {
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
      const double threshold = predictor.GetKs(config, jump) * std::pow(predictor.GetD(config, jump), 2) / 12;
      const double model_barrier = predictor.GetModelBarrierFromLatticeIdPair(config, jump);
      const auto [barrier, energy_change] = predictor.GetBarrierAndDiffFromLatticeIdPair(config, jump);
      ASSERT_EQ(energy_change, dE);
      ASSERT_EQ(model_barrier, GetQuarticModelBarrier(dE, predictor.GetD(config, jump), predictor.GetKs(config, jump)));
      // Leave out a jump that sits on the threshold to within the resolution of the comparison with dE.
      if (std::abs(std::abs(dE) - threshold) < 1e-4 * threshold) {
        continue;
      }
      if (dE < -threshold) {
        ++num_without_forward_barrier;
        ASSERT_LT(model_barrier, 0.0) << "step " << step;
        ASSERT_EQ(barrier, 0.0) << "step " << step;
      } else if (dE > threshold) {
        ++num_without_backward_barrier;
        ASSERT_LT(model_barrier, dE) << "step " << step;
        ASSERT_EQ(barrier, dE) << "step " << step;
      } else {
        ++num_regular;
        ASSERT_GT(model_barrier, std::max(dE, 0.0)) << "step " << step;
        ASSERT_EQ(barrier, model_barrier) << "step " << step;
      }
    }
    config.LatticeJump({vacancy, config.GetFirstNeighborsAdjacencyList()[vacancy][pick(generator)]});
  }
  // The model was made so that all three cases occur.
  EXPECT_GE(num_regular, 50u);
  EXPECT_GE(num_without_forward_barrier, 50u);
  EXPECT_GE(num_without_backward_barrier, 50u);
}

// Jumping back passes the same saddle point, or none: the two barriers differ by the energy change whether the model
// has a barrier or not.
TEST_F(JumpsWithoutBarrier, ReverseJumpIsConsistentWithAndWithoutABarrier) {
  const VacancyMigrationPredictorQuartic predictor(model_filename_, config_, TestElementSet());
  auto config = config_;
  std::mt19937_64 generator(11);
  std::uniform_int_distribution<size_t> pick(0, constants::kNumFirstNearestNeighbors - 1);
  size_t num_limited = 0;
  for (size_t step = 0; step < 200; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    const size_t target = config.GetFirstNeighborsAdjacencyList()[vacancy][pick(generator)];
    const auto [forward_barrier, forward_change] =
        predictor.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, target});
    ASSERT_GE(forward_barrier, std::max(forward_change, 0.0)) << "step " << step;
    config.LatticeJump({vacancy, target});
    const auto [backward_barrier, backward_change] =
        predictor.GetBarrierAndDiffFromLatticeIdPair(config, {target, vacancy});
    ASSERT_NEAR(backward_change, -forward_change, 1e-12) << "step " << step;
    ASSERT_NEAR(backward_barrier, forward_barrier - forward_change, 1e-9) << "step " << step;
    ASSERT_GE(backward_barrier, 0.0) << "step " << step;
    if (forward_barrier == 0.0 || forward_barrier == forward_change) {
      ++num_limited;
      // exactly one of the two directions is free
      ASSERT_NEAR(std::min(forward_barrier, backward_barrier), 0.0, 1e-9) << "step " << step;
      ASSERT_NEAR(std::max(forward_barrier, backward_barrier), std::abs(forward_change), 1e-9) << "step " << step;
    }
  }
  EXPECT_GE(num_limited, 20u);
}

TEST_F(JumpsWithoutBarrier, CachedPredictorAgreesWithThePlainPredictor) {
  const VacancyMigrationPredictorQuartic plain(model_filename_, config_, TestElementSet());
  const VacancyMigrationPredictorQuarticLru cached(model_filename_, config_, TestElementSet(), 64);
  auto config = config_;
  size_t num_limited = 0;
  for (size_t step = 0; step < 20; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    for (const size_t neighbor: config.GetFirstNeighborsAdjacencyList()[vacancy]) {
      const auto expected = plain.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, neighbor});
      if (expected.first == std::max(expected.second, 0.0)) {
        ++num_limited;
      }
      // The first call fills the cache, the second one is answered from it.
      for (int repetition = 0; repetition < 2; ++repetition) {
        const auto actual = cached.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, neighbor});
        ASSERT_EQ(actual.first, expected.first) << "step " << step;
        ASSERT_EQ(actual.second, expected.second) << "step " << step;
      }
      ASSERT_EQ(cached.GetModelBarrierFromLatticeIdPair(config, {vacancy, neighbor}),
                plain.GetModelBarrierFromLatticeIdPair(config, {vacancy, neighbor}));
    }
    config.LatticeJump({vacancy, config.GetFirstNeighborsAdjacencyList()[vacancy][step % 12]});
  }
  EXPECT_GE(num_limited, 50u);
}
}    // namespace pred

#include "EnergyChangePredictorPairSite.h"
#include "EnergyPredictor.h"
#include "VacancyMigrationPredictorQuartic.h"
#include "VacancyMigrationPredictorQuarticLru.h"
#include "cfg/test/TemporaryPath.hpp"
#include "pred/test/SyntheticModel.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <random>
#include <stdexcept>

namespace pred {
// The Monte Carlo methods never evaluate the total energy while they run: KMC and CMC add up the energy changes of
// the incremental predictors. These tests tie the incremental predictors to EnergyPredictor and to each other.
class PredictorConsistency : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    temporary_directory_ = std::make_unique<test::TemporaryDirectory>();
    directory_ = temporary_directory_->GetPath();
    config_ = MakeTestAlloy();
    model_ = MakeSyntheticModel(config_, TestElementSet());
    model_filename_ = WriteModel(directory_ / "model.json", model_);
  }

  static void TearDownTestSuite() {
    temporary_directory_.reset();
  }

  static inline std::unique_ptr<test::TemporaryDirectory> temporary_directory_{};
  static inline std::filesystem::path directory_{};
  static inline cfg::Config config_{};
  static inline nlohmann::json model_{};
  static inline std::string model_filename_{};
};

TEST_F(PredictorConsistency, SwapEnergyChangeEqualsTotalEnergyDifference) {
  const EnergyPredictor total_energy_predictor(model_filename_, TestElementSet());
  const EnergyChangePredictorPairSite energy_change_predictor(model_filename_, config_, TestElementSet());
  auto config = config_;
  std::mt19937_64 generator(3);
  std::uniform_int_distribution<size_t> pick(0, config.GetNumAtoms() - 1);
  double energy = total_energy_predictor.GetEnergy(config);
  size_t num_neighboring = 0;
  size_t num_distant = 0;
  for (size_t trial = 0; trial < 400 && num_neighboring + num_distant < 60; ++trial) {
    // Sites that are neighbours share clusters and take a different code path than distant sites.
    const size_t first = pick(generator);
    size_t second = pick(generator);
    if (trial % 4 == 1) {
      second = config.GetFirstNeighborsAdjacencyList()[first][pick(generator) % constants::kNumFirstNearestNeighbors];
    } else if (trial % 4 == 2) {
      second = config.GetSecondNeighborsAdjacencyList()[first][pick(generator) % constants::kNumSecondNearestNeighbors];
    } else if (trial % 4 == 3) {
      second = config.GetThirdNeighborsAdjacencyList()[first][pick(generator) % constants::kNumThirdNearestNeighbors];
    }
    if (config.GetElementAtLatticeId(first) == config.GetElementAtLatticeId(second)) {
      continue;
    }
    (config.FindDistanceLabelBetweenLattice(first, second) > 0 ? num_neighboring : num_distant)++;
    const double predicted_change = energy_change_predictor.GetDeFromLatticeIdPair(config, {first, second});
    config.LatticeJump({first, second});
    const double new_energy = total_energy_predictor.GetEnergy(config);
    ASSERT_NEAR(predicted_change, new_energy - energy, 1e-9) << "sites " << first << " and " << second;
    energy = new_energy;
  }
  EXPECT_GE(num_neighboring, 10u);
  EXPECT_GE(num_distant, 10u);
}

// The model has no terms for two vacancies in one cluster. Whichever of the two sites is named first, a swap of the
// vacancy with a neighbour must not be evaluated through a state in which both sites are vacant.
TEST_F(PredictorConsistency, SwapOfTheVacancyWithANeighbourWorksInBothOrders) {
  const EnergyPredictor total_energy_predictor(model_filename_, TestElementSet());
  const EnergyChangePredictorPairSite energy_change_predictor(model_filename_, config_, TestElementSet());
  const size_t vacancy = config_.GetVacancyLatticeId();
  const double energy = total_energy_predictor.GetEnergy(config_);
  size_t num_neighbors = 0;
  for (const auto *adjacency_list: {&config_.GetFirstNeighborsAdjacencyList(),
                                    &config_.GetSecondNeighborsAdjacencyList(),
                                    &config_.GetThirdNeighborsAdjacencyList()}) {
    for (const auto neighbor: (*adjacency_list)[vacancy]) {
      auto swapped_config = config_;
      swapped_config.LatticeJump({vacancy, neighbor});
      const double expected_change = total_energy_predictor.GetEnergy(swapped_config) - energy;
      EXPECT_NEAR(energy_change_predictor.GetDeFromLatticeIdPair(config_, {vacancy, neighbor}), expected_change, 1e-9);
      EXPECT_NEAR(energy_change_predictor.GetDeFromLatticeIdPair(config_, {neighbor, vacancy}), expected_change, 1e-9);
      ++num_neighbors;
    }
  }
  EXPECT_EQ(num_neighbors, constants::kNumFirstNearestNeighbors + constants::kNumSecondNearestNeighbors +
      constants::kNumThirdNearestNeighbors);
}

TEST_F(PredictorConsistency, JumpEnergyChangeEqualsTotalEnergyDifferenceAndReverseJumpIsConsistent) {
  const EnergyPredictor total_energy_predictor(model_filename_, TestElementSet());
  const VacancyMigrationPredictorQuartic barrier_predictor(model_filename_, config_, TestElementSet());
  auto config = config_;
  std::mt19937_64 generator(5);
  std::uniform_int_distribution<size_t> pick(0, constants::kNumFirstNearestNeighbors - 1);
  double energy = total_energy_predictor.GetEnergy(config);
  for (size_t step = 0; step < 60; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    const size_t target = config.GetFirstNeighborsAdjacencyList()[vacancy][pick(generator)];
    const auto [forward_barrier, forward_change] =
        barrier_predictor.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, target});
    config.LatticeJump({vacancy, target});
    const double new_energy = total_energy_predictor.GetEnergy(config);
    ASSERT_NEAR(forward_change, new_energy - energy, 1e-9) << "step " << step;

    // Jumping back has to undo the energy change and to pass the same saddle point.
    const auto [backward_barrier, backward_change] =
        barrier_predictor.GetBarrierAndDiffFromLatticeIdPair(config, {target, vacancy});
    ASSERT_NEAR(backward_change, -forward_change, 1e-12) << "step " << step;
    ASSERT_NEAR(backward_barrier, forward_barrier - forward_change, 1e-9) << "step " << step;
    ASSERT_TRUE(std::isfinite(forward_barrier));
    energy = new_energy;
  }
}

TEST_F(PredictorConsistency, CachedPredictorReturnsTheBarriersOfThePlainPredictor) {
  const VacancyMigrationPredictorQuartic plain(model_filename_, config_, TestElementSet());
  const VacancyMigrationPredictorQuarticLru cached(model_filename_, config_, TestElementSet(), 64);
  auto config = config_;
  for (size_t step = 0; step < 20; ++step) {
    const size_t vacancy = config.GetVacancyLatticeId();
    for (const size_t neighbor: config.GetFirstNeighborsAdjacencyList()[vacancy]) {
      const auto expected = plain.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, neighbor});
      // The first call fills the cache, the second one is answered from it.
      for (int repetition = 0; repetition < 2; ++repetition) {
        const auto actual = cached.GetBarrierAndDiffFromLatticeIdPair(config, {vacancy, neighbor});
        ASSERT_DOUBLE_EQ(actual.first, expected.first);
        ASSERT_DOUBLE_EQ(actual.second, expected.second);
      }
    }
    config.LatticeJump({vacancy, config.GetFirstNeighborsAdjacencyList()[vacancy][step % 12]});
  }
}

namespace {
// Gives a test access to the two fitted quantities of the barrier model.
class ExposedBarrierPredictor : public VacancyMigrationPredictorQuartic {
 public:
  using VacancyMigrationPredictorQuartic::VacancyMigrationPredictorQuartic;
  using VacancyMigrationPredictorQuartic::GetD;
  using VacancyMigrationPredictorQuartic::GetKs;
};
}    // namespace

// The barrier of a jump is a function of the elements around it. Two jumps elsewhere in the lattice, in any of the
// twelve directions, with the same elements at the same places of their neighbourhoods have to get the same barrier
// to the last bit: the cache hands the barrier of one to the other.
TEST_F(PredictorConsistency, JumpsWithTheSameNeighbourhoodHaveTheSameBarrierToTheLastBit) {
  // Without the offset and the scale of the fitted quantities, which would round a difference in the last bit of the
  // sum over the clusters away in most cases.
  auto model = model_;
  for (const auto &element: TestElementSet()) {
    auto &parameters = model[element.GetString()];
    parameters["mu_D"] = 0.0;
    parameters["mu_Ks"] = 0.0;
    parameters["sigma_D"] = 1.0;
    parameters["sigma_Ks"] = 1.0;
  }
  const auto model_filename = WriteModel(directory_ / "model_without_output_scale.json", model);
  const auto pure_config = cfg::GenerateFCC({8, 8, 8}, Element("Al"));
  const ExposedBarrierPredictor barrier_predictor(model_filename, pure_config, TestElementSet());
  const std::vector<Element> elements(TestElementSet().begin(), TestElementSet().end());
  // The elements of the neighbourhood in the order of the state list of a pair, with the vacancy on the first site.
  const auto decorate = [&](const std::pair<size_t, size_t> &pair, const std::vector<Element> &neighbourhood) {
    auto config = pure_config;
    const auto state_list = GetSortedLatticeVectorStateOfPair(pure_config, pair);
    for (size_t k = 0; k < state_list.size(); ++k) {
      const auto lattice_id = state_list[k].GetId();
      config.SetAtomElementTypeAtLattice(lattice_id, lattice_id == pair.first ? Element("X") : neighbourhood[k]);
    }
    return config;
  };
  std::mt19937_64 generator(17);
  const std::pair<size_t, size_t> reference_pair{0, pure_config.GetFirstNeighborsAdjacencyList()[0][0]};
  size_t num_compared = 0;
  for (int trial = 0; trial < 5; ++trial) {
    std::vector<Element> neighbourhood;
    for (size_t k = 0; k < constants::kNumThirdNearestSetSizeOfPair; ++k) {
      neighbourhood.push_back(elements[generator() % elements.size()]);
    }
    const auto reference_config = decorate(reference_pair, neighbourhood);
    const auto expected = barrier_predictor.GetBarrierAndDiffFromLatticeIdPair(reference_config, reference_pair);
    const double expected_d = barrier_predictor.GetD(reference_config, reference_pair);
    const double expected_ks = barrier_predictor.GetKs(reference_config, reference_pair);
    ASSERT_TRUE(std::isfinite(expected.first));
    for (const size_t first: {size_t{0}, size_t{777}, size_t{1500}}) {
      for (const size_t second: pure_config.GetFirstNeighborsAdjacencyList()[first]) {
        const std::pair<size_t, size_t> pair{first, second};
        const auto config = decorate(pair, neighbourhood);
        const auto actual = barrier_predictor.GetBarrierAndDiffFromLatticeIdPair(config, pair);
        ASSERT_EQ(barrier_predictor.GetD(config, pair), expected_d) << "trial " << trial << ", jump " << first << " -> "
                                                                    << second;
        ASSERT_EQ(barrier_predictor.GetKs(config, pair), expected_ks) << "trial " << trial << ", jump " << first
                                                                      << " -> " << second;
        ASSERT_EQ(actual.first, expected.first) << "trial " << trial << ", jump " << first << " -> " << second;
        ASSERT_EQ(actual.second, expected.second) << "trial " << trial << ", jump " << first << " -> " << second;
        ++num_compared;
      }
    }
  }
  EXPECT_EQ(num_compared, 5u * 3u * constants::kNumFirstNearestNeighbors);
}

// A box that is shorter than three cutoff lengths used to get neighbour lists with repeated entries.
TEST_F(PredictorConsistency, EnergyPerAtomOfThePureSolventDoesNotDependOnTheBoxSize) {
  const EnergyPredictor total_energy_predictor(model_filename_, TestElementSet());
  const auto energy_per_atom = [&](const size_t factor) {
    const auto config = cfg::GenerateFCC({factor, factor, factor}, Element("Al"));
    return total_energy_predictor.GetEnergy(config) / static_cast<double>(config.GetNumAtoms());
  };
  const double reference = energy_per_atom(5);
  EXPECT_NEAR(energy_per_atom(3), reference, 1e-10);
  EXPECT_NEAR(energy_per_atom(4), reference, 1e-10);
}

// Release builds switch the size checks of Eigen off, so a file for other elements gave meaningless energies.
TEST_F(PredictorConsistency, CoefficientsForAnotherElementSetAreRejected) {
  const std::set<Element> binary{Element("Al"), Element("Mg")};
  EXPECT_THROW(EnergyPredictor(model_filename_, binary), std::runtime_error);
  EXPECT_THROW(EnergyChangePredictorPairSite(model_filename_, config_, binary), std::runtime_error);
  EXPECT_THROW(VacancyMigrationPredictorQuartic(model_filename_, config_, binary), std::runtime_error);
}

TEST_F(PredictorConsistency, InconsistentCoefficientFilesAreRejected) {
  auto truncated = model_;
  truncated["Al"]["mu_x_mmm"].erase(0);
  EXPECT_THROW(VacancyMigrationPredictorQuartic(WriteModel(directory_ / "truncated.json", truncated), config_,
                                                TestElementSet()),
               std::runtime_error);

  auto longer = model_;
  longer["Mg"]["theta_Ks"].push_back(0.0);
  EXPECT_THROW(
      VacancyMigrationPredictorQuartic(WriteModel(directory_ / "longer.json", longer), config_, TestElementSet()),
      std::runtime_error);

  auto incomplete = model_;
  incomplete.erase("Zn");
  EXPECT_THROW(VacancyMigrationPredictorQuartic(WriteModel(directory_ / "incomplete.json", incomplete), config_,
                                                TestElementSet()),
               std::runtime_error);

  auto no_base = model_;
  no_base.erase("Base");
  EXPECT_THROW(EnergyPredictor(WriteModel(directory_ / "no_base.json", no_base), TestElementSet()),
               std::runtime_error);
}
}    // namespace pred

#include "CanonicalMcSerial.h"
#include "EnergyPredictor.h"
#include "KineticMcFirstSerial.h"
#include "KineticMcSecondSerial.h"
#include "VacancyMigrationPredictorQuartic.h"
#include "cfg/test/TemporaryPath.hpp"
#include "pred/test/SyntheticModel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mc {
namespace {
// The simulations write their log and configurations to the working directory.
class WorkingDirectory {
 public:
  explicit WorkingDirectory(const std::filesystem::path &directory)
      : previous_(std::filesystem::current_path()) {
    std::filesystem::create_directories(directory);
    std::filesystem::current_path(directory);
  }
  WorkingDirectory(const WorkingDirectory &) = delete;
  WorkingDirectory &operator=(const WorkingDirectory &) = delete;
  ~WorkingDirectory() {
    std::filesystem::current_path(previous_);
  }

 private:
  const std::filesystem::path previous_;
};

std::string ReadFile(const std::filesystem::path &filename) {
  std::ifstream ifs(filename);
  std::stringstream buffer;
  buffer << ifs.rdbuf();
  return buffer.str();
}

// Collects what is written to std::cerr while it exists.
class CapturedErrorOutput {
 public:
  CapturedErrorOutput()
      : previous_(std::cerr.rdbuf(buffer_.rdbuf())) {}
  CapturedErrorOutput(const CapturedErrorOutput &) = delete;
  CapturedErrorOutput &operator=(const CapturedErrorOutput &) = delete;
  ~CapturedErrorOutput() {
    std::cerr.rdbuf(previous_);
  }
  [[nodiscard]] std::string GetString() const {
    return buffer_.str();
  }

 private:
  std::ostringstream buffer_{};
  std::streambuf *const previous_;
};

// The rows of a log file without its header, split into columns.
std::vector<std::vector<std::string>> ReadLog(const std::filesystem::path &filename) {
  std::ifstream ifs(filename);
  std::vector<std::vector<std::string>> rows;
  std::string line;
  std::getline(ifs, line);
  while (std::getline(ifs, line)) {
    std::istringstream line_stream(line);
    std::vector<std::string> columns;
    for (std::string column; std::getline(line_stream, column, '\t');) {
      columns.push_back(column);
    }
    rows.push_back(columns);
  }
  return rows;
}

// Lets a test look at single steps of KineticMcSecondSerial.
class SecondOrderStepper : public KineticMcSecondSerial {
 public:
  using KineticMcSecondSerial::KineticMcSecondSerial;
  // Evaluates the next step without making it and returns the time it takes.
  double EvaluateStep() {
    BuildEventList();
    return CalculateTime();
  }
  // The probabilities that the evaluated step ends on each of the twelve neighbor sites of the vacancy.
  [[nodiscard]] std::array<double, kEventListSize> GetStepProbabilities() const {
    std::array<double, kEventListSize> probabilities{};
    for (size_t i = 0; i < kEventListSize; ++i) {
      probabilities[i] = event_k_i_list_[i].GetProbability();
    }
    return probabilities;
  }
  void MakeStep() {
    OneStepSimulation();
  }
  [[nodiscard]] const cfg::Config &GetConfig() const {
    return config_;
  }
  [[nodiscard]] size_t GetVacancyLatticeId() const {
    return vacancy_lattice_id_;
  }
  // The site the vacancy came from.
  [[nodiscard]] size_t GetPreviousLatticeId() const {
    return previous_j_lattice_id_;
  }
};
}    // namespace

class Simulation : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    temporary_directory_ = std::make_unique<test::TemporaryDirectory>();
    directory_ = temporary_directory_->GetPath();
    config_ = pred::MakeTestAlloy();
    model_filename_ =
        pred::WriteModel(directory_ / "model.json", pred::MakeSyntheticModel(config_, pred::TestElementSet()));
    omp_set_num_threads(2);
  }

  static void TearDownTestSuite() {
    temporary_directory_.reset();
  }

  // 300 steps of first-order KMC at 600 K, logged every 10 steps.
  static void RunKineticMc(const std::string &name, const unsigned long long int seed) {
    const WorkingDirectory working_directory(directory_ / name);
    KineticMcFirstSerial(config_, 10, 100, 300, 0, 0, 0, 0, 600, pred::TestElementSet(), model_filename_, "", false,
                         true, seed)
        .Simulate();
  }

  // 300 steps of second-order KMC at 600 K, logged every 10 steps.
  static void RunKineticMcSecond(const std::string &name, const unsigned long long int seed) {
    const WorkingDirectory working_directory(directory_ / name);
    KineticMcSecondSerial(config_, 10, 100, 300, 0, 0, 0, 0, 600, pred::TestElementSet(), model_filename_, "", false,
                          true, seed)
        .Simulate();
  }

  // 2000 trial swaps at 600 K with every step logged and an averaging window of 50 steps.
  static void RunCanonicalMc(const std::string &name, const unsigned long long int seed) {
    const WorkingDirectory working_directory(directory_ / name);
    CanonicalMcSerial(config_, 1, 1000, 2000, 50, 0, 0, 600, pred::TestElementSet(), model_filename_, seed)
        .Simulate();
  }

  static inline std::unique_ptr<test::TemporaryDirectory> temporary_directory_{};
  static inline std::filesystem::path directory_{};
  static inline cfg::Config config_{};
  static inline std::string model_filename_{};
};

TEST_F(Simulation, KineticMcIsReproducibleForAGivenSeed) {
  RunKineticMc("kmc_seed_11", 11);
  RunKineticMc("kmc_seed_11_again", 11);
  RunKineticMc("kmc_seed_12", 12);
  const auto log = ReadFile(directory_ / "kmc_seed_11" / "kmc_log.txt");
  EXPECT_GT(ReadLog(directory_ / "kmc_seed_11" / "kmc_log.txt").size(), 30u);
  EXPECT_EQ(log, ReadFile(directory_ / "kmc_seed_11_again" / "kmc_log.txt"));
  EXPECT_NE(log, ReadFile(directory_ / "kmc_seed_12" / "kmc_log.txt"));
}

// The energy column is the sum of the predicted energy changes of all jumps.
TEST_F(Simulation, KineticMcEnergyColumnFollowsTheTotalEnergy) {
  RunKineticMc("kmc_energy", 21);
  const auto rows = ReadLog(directory_ / "kmc_energy" / "kmc_log.txt");
  ASSERT_FALSE(rows.empty());
  ASSERT_EQ(rows.back().at(0), "300");
  const pred::EnergyPredictor total_energy_predictor(model_filename_, pred::TestElementSet());
  const auto final_config = cfg::Config::ReadConfig((directory_ / "kmc_energy" / "end.cfg.gz").string());
  EXPECT_NEAR(std::stod(rows.back().at(3)),
              total_energy_predictor.GetEnergy(final_config) - total_energy_predictor.GetEnergy(config_), 1e-8);
  // time and energy of every row are finite numbers
  for (const auto &row: rows) {
    ASSERT_TRUE(std::isfinite(std::stod(row.at(1)))) << "step " << row.at(0);
    ASSERT_TRUE(std::isfinite(std::stod(row.at(3)))) << "step " << row.at(0);
  }
}

TEST_F(Simulation, KineticMcSecondSerialIsReproducibleForAGivenSeed) {
  RunKineticMcSecond("kmc_second_seed_11", 11);
  RunKineticMcSecond("kmc_second_seed_11_again", 11);
  RunKineticMcSecond("kmc_second_seed_12", 12);
  const auto log = ReadFile(directory_ / "kmc_second_seed_11" / "kmc_log.txt");
  EXPECT_GT(ReadLog(directory_ / "kmc_second_seed_11" / "kmc_log.txt").size(), 30u);
  EXPECT_EQ(log, ReadFile(directory_ / "kmc_second_seed_11_again" / "kmc_log.txt"));
  EXPECT_NE(log, ReadFile(directory_ / "kmc_second_seed_12" / "kmc_log.txt"));
}

TEST_F(Simulation, KineticMcSecondSerialEnergyColumnFollowsTheTotalEnergy) {
  RunKineticMcSecond("kmc_second_energy", 21);
  const auto rows = ReadLog(directory_ / "kmc_second_energy" / "kmc_log.txt");
  ASSERT_FALSE(rows.empty());
  ASSERT_EQ(rows.back().at(0), "300");
  const pred::EnergyPredictor total_energy_predictor(model_filename_, pred::TestElementSet());
  const auto final_config = cfg::Config::ReadConfig((directory_ / "kmc_second_energy" / "end.cfg.gz").string());
  EXPECT_NEAR(std::stod(rows.back().at(3)),
              total_energy_predictor.GetEnergy(final_config) - total_energy_predictor.GetEnergy(config_), 1e-8);
  // The time never decreases, and time and energy of every row are finite numbers.
  double previous_time = 0.0;
  for (const auto &row: rows) {
    const double time = std::stod(row.at(1));
    ASSERT_TRUE(std::isfinite(time)) << "step " << row.at(0);
    ASSERT_TRUE(std::isfinite(std::stod(row.at(3)))) << "step " << row.at(0);
    ASSERT_GE(time, previous_time) << "step " << row.at(0);
    previous_time = time;
  }
}

// In a second-order step the vacancy on site k jumps to a neighbor site i and back as often as it happens to, and the
// step ends on the first site i from which it does not return at once. Only the first jump of the step cannot lead to
// j, the site the vacancy came from. The method has closed formulas for the probability that the step ends on a site
// and for its mean time. This test follows the same process attempt by attempt, with the barriers taken straight from
// the model, and compares.
TEST_F(Simulation, KineticMcSecondSerialStepAgreesWithTheSumOverItsPaths) {
  const WorkingDirectory working_directory(directory_ / "kmc_second_paths");
  const double temperature = 600;
  const double beta = 1.0 / constants::kBoltzmann / temperature;
  const size_t num_neighbors = constants::kNumFirstNearestNeighbors;
  SecondOrderStepper kinetic_mc(config_, 1, 100, 100, 0, 0, 0, 0, temperature, pred::TestElementSet(), model_filename_,
                                "", false, false, 61);
  const pred::VacancyMigrationPredictorQuartic predictor(model_filename_, config_, pred::TestElementSet());
  const auto get_rate = [&](const cfg::Config &config, const size_t from, const size_t to) {
    return std::exp(-predictor.GetBarrierAndDiffFromLatticeIdPair(config, {from, to}).first * beta);
  };
  for (size_t step = 0; step < 40; ++step) {
    const double step_time = kinetic_mc.EvaluateStep();
    const auto step_probabilities = kinetic_mc.GetStepProbabilities();

    auto config = kinetic_mc.GetConfig();
    const size_t k = kinetic_mc.GetVacancyLatticeId();
    const size_t j = kinetic_mc.GetPreviousLatticeId();
    const auto neighbors_of_k = config.GetFirstNeighborsAdjacencyList()[k];
    // jump_probability: of the jump k -> i among the twelve jumps of the vacancy on k,
    // return_probability: of the jump i -> k among the twelve jumps of the vacancy on i,
    // residence_time: the mean time the vacancy stays on i.
    std::vector<double> jump_probability(num_neighbors), return_probability(num_neighbors);
    std::vector<double> residence_time(num_neighbors);
    double total_rate_k = 0.0;
    double jump_probability_to_j = 0.0;
    for (size_t i = 0; i < num_neighbors; ++i) {
      jump_probability[i] = get_rate(config, k, neighbors_of_k[i]);
      total_rate_k += jump_probability[i];
    }
    for (size_t i = 0; i < num_neighbors; ++i) {
      const size_t site_i = neighbors_of_k[i];
      jump_probability[i] /= total_rate_k;
      if (site_i == j) {
        jump_probability_to_j = jump_probability[i];
      }
      config.LatticeJump({k, site_i});
      double total_rate_i = 0.0;
      for (const auto l: config.GetFirstNeighborsAdjacencyList()[site_i]) {
        total_rate_i += get_rate(config, site_i, l);
      }
      return_probability[i] = get_rate(config, site_i, k) / total_rate_i;
      residence_time[i] = 1.0 / total_rate_i / constants::kPrefactor;
      config.LatticeJump({site_i, k});
    }
    const double residence_time_k = 1.0 / total_rate_k / constants::kPrefactor;

    // `weight` is the probability that the vacancy is on k for the given attempt, that is, that it has returned from
    // every jump of the step so far.
    std::vector<double> expected_probabilities(num_neighbors, 0.0);
    double expected_time = 0.0;
    double weight = 1.0;
    for (size_t attempt = 0; attempt < 100000 && weight > 1e-18; ++attempt) {
      expected_time += weight * residence_time_k;
      double weight_of_returns = 0.0;
      for (size_t i = 0; i < num_neighbors; ++i) {
        double probability = jump_probability[i];
        if (attempt == 0) {
          probability = neighbors_of_k[i] == j ? 0.0 : jump_probability[i] / (1.0 - jump_probability_to_j);
        }
        expected_probabilities[i] += weight * probability * (1.0 - return_probability[i]);
        expected_time += weight * probability * return_probability[i] * residence_time[i];
        weight_of_returns += weight * probability * return_probability[i];
      }
      weight = weight_of_returns;
    }
    ASSERT_LT(weight, 1e-12) << "step " << step;

    EXPECT_NEAR(step_time, expected_time, 1e-9 * expected_time) << "step " << step;
    double sum_of_probabilities = 0.0;
    for (size_t i = 0; i < num_neighbors; ++i) {
      EXPECT_NEAR(step_probabilities[i], expected_probabilities[i], 1e-9) << "step " << step << ", neighbor " << i;
      sum_of_probabilities += step_probabilities[i];
    }
    EXPECT_NEAR(sum_of_probabilities, 1.0, 1e-12) << "step " << step;
    // A step takes at least as long as the first jump away from k.
    EXPECT_GE(step_time, residence_time_k * (1.0 - 1e-12)) << "step " << step;

    kinetic_mc.MakeStep();
  }
}

TEST_F(Simulation, CanonicalMcIsReproducibleForAGivenSeed) {
  RunCanonicalMc("cmc_seed_11", 11);
  RunCanonicalMc("cmc_seed_11_again", 11);
  RunCanonicalMc("cmc_seed_12", 12);
  const auto log = ReadFile(directory_ / "cmc_seed_11" / "cmc_log.txt");
  EXPECT_EQ(log, ReadFile(directory_ / "cmc_seed_11_again" / "cmc_log.txt"));
  EXPECT_NE(log, ReadFile(directory_ / "cmc_seed_12" / "cmc_log.txt"));
}

TEST_F(Simulation, CanonicalMcLogsEnergyAndItsPlainAverage) {
  RunCanonicalMc("cmc_energy", 31);
  const auto rows = ReadLog(directory_ / "cmc_energy" / "cmc_log.txt");
  ASSERT_EQ(rows.size(), 2001u);    // steps 0 to 2000
  // columns: steps, temperature, energy, average_energy, absolute_energy
  const size_t window = 50;
  double sum = 0;
  for (size_t step = 0; step < rows.size(); ++step) {
    sum += std::stod(rows[step].at(2));
    if (step >= window) {
      sum -= std::stod(rows[step - window].at(2));
    }
    ASSERT_NEAR(std::stod(rows[step].at(3)), sum / static_cast<double>(std::min(step + 1, window)), 1e-9)
        << "step " << step;
  }
  const pred::EnergyPredictor total_energy_predictor(model_filename_, pred::TestElementSet());
  const auto final_config = cfg::Config::ReadConfig((directory_ / "cmc_energy" / "end.cfg.gz").string());
  EXPECT_NEAR(std::stod(rows.back().at(2)),
              total_energy_predictor.GetEnergy(final_config) - total_energy_predictor.GetEnergy(config_), 1e-8);
  EXPECT_NEAR(std::stod(rows.back().at(4)), total_energy_predictor.GetEnergy(final_config), 1e-8);
}

// On a model whose energy changes are so large that many jumps have a negative barrier in one direction, KMC has to
// take the barriers as they are and to tell how many of the jumps it took had a negative one.
TEST_F(Simulation, KineticMcCountsTheNegativeBarriers) {
  const auto model_filename =
      pred::WriteModel(directory_ / "model_with_large_energy_changes.json",
                       pred::MakeSyntheticModelWithScaledEnergies(config_, pred::TestElementSet(), 300));
  std::string warnings;
  {
    const WorkingDirectory working_directory(directory_ / "kmc_negative_barriers");
    const CapturedErrorOutput captured;
    // 1000 steps, each of them logged, at a temperature at which the vacancy climbs about as often as it descends.
    KineticMcFirstSerial(config_, 1, 100, 1000, 0, 0, 0, 0, 20000, pred::TestElementSet(), model_filename, "", false,
                         false, 41)
        .Simulate();
    warnings = captured.GetString();
  }
  const auto rows = ReadLog(directory_ / "kmc_negative_barriers" / "kmc_log.txt");
  ASSERT_EQ(rows.size(), 1001u);    // steps 0 to 1000
  // columns: steps, time, temperature, energy, Ea, dE, ...
  size_t num_negative_forward = 0;
  size_t num_negative_backward = 0;
  std::string first_step;
  for (const auto &row: rows) {
    const double barrier = std::stod(row.at(4));
    ASSERT_TRUE(std::isfinite(std::stod(row.at(1)))) << "step " << row.at(0);
    if (barrier < 0.0) {
      ++num_negative_forward;
    } else if (barrier < std::stod(row.at(5))) {
      ++num_negative_backward;
    } else {
      continue;
    }
    if (first_step.empty()) {
      first_step = row.at(0);
    }
  }
  EXPECT_GE(num_negative_forward, 5u);
  EXPECT_GE(num_negative_backward, 5u);
  // The first one is reported when it happens, because a long run is stopped from outside and never gets to the
  // summary at its end.
  ASSERT_FALSE(first_step.empty());
  const std::string first_warning = "Warning: step " + first_step + ": so far 1 of " +
      std::to_string(std::stoull(first_step) + 1) + " jumps had a negative barrier";
  EXPECT_EQ(warnings.rfind(first_warning, 0), 0u) << warnings;
  const std::string summary = "Warning: in this run " + std::to_string(num_negative_forward + num_negative_backward) +
      " of 1001 jumps had a negative barrier (" + std::to_string(num_negative_forward) + " forward, " +
      std::to_string(num_negative_backward) + " only backward)";
  EXPECT_NE(warnings.find(summary), std::string::npos) << warnings;
}

TEST_F(Simulation, KineticMcWritesNoWarningWhileEveryBarrierIsPositive) {
  std::string warnings;
  {
    const CapturedErrorOutput captured;
    RunKineticMc("kmc_positive_barriers", 51);
    warnings = captured.GetString();
  }
  EXPECT_EQ(warnings, "");
  for (const auto &row: ReadLog(directory_ / "kmc_positive_barriers" / "kmc_log.txt")) {
    const double barrier = std::stod(row.at(4));
    const double energy_change = std::stod(row.at(5));
    ASSERT_GT(barrier, std::max(energy_change, 0.0)) << "step " << row.at(0);
  }
}

// KMC follows one vacancy and its barrier model has no terms for a second one.
TEST_F(Simulation, KineticMcNeedsExactlyOneVacancy) {
  const WorkingDirectory working_directory(directory_ / "kmc_vacancies");
  auto no_vacancy = config_;
  no_vacancy.SetAtomElementTypeAtLattice(0, Element("Al"));
  auto two_vacancies = config_;
  two_vacancies.SetAtomElementTypeAtLattice(100, Element("X"));
  EXPECT_THROW(KineticMcFirstSerial(no_vacancy, 10, 100, 300, 0, 0, 0, 0, 600, pred::TestElementSet(), model_filename_,
                                    "", false, false, 1),
               std::runtime_error);
  EXPECT_THROW(KineticMcFirstSerial(two_vacancies, 10, 100, 300, 0, 0, 0, 0, 600, pred::TestElementSet(),
                                    model_filename_, "", false, false, 1),
               std::runtime_error);
  EXPECT_THROW(KineticMcSecondSerial(two_vacancies, 10, 100, 300, 0, 0, 0, 0, 600, pred::TestElementSet(),
                                     model_filename_, "", false, false, 1),
               std::runtime_error);
}

TEST_F(Simulation, ElementsOutsideElementSetAreRejected) {
  const WorkingDirectory working_directory(directory_ / "foreign_element");
  auto with_copper = config_;
  with_copper.SetAtomElementTypeAtLattice(100, Element("Cu"));
  EXPECT_THROW(
      CanonicalMcSerial(with_copper, 1, 1000, 2000, 50, 0, 0, 600, pred::TestElementSet(), model_filename_, 1),
      std::runtime_error);
  EXPECT_THROW(KineticMcFirstSerial(with_copper, 10, 100, 300, 0, 0, 0, 0, 600, pred::TestElementSet(),
                                    model_filename_, "", false, false, 1),
               std::runtime_error);
}
}    // namespace mc

#include "CanonicalMcSerial.h"
#include "EnergyPredictor.h"
#include "KineticMcFirstSerial.h"
#include "cfg/test/TemporaryPath.hpp"
#include "pred/test/SyntheticModel.hpp"

#include <algorithm>
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

// On a model whose energy changes are so large that four of ten jumps have no barrier in the quartic form, KMC has
// to run without a negative barrier and to tell how many of the jumps it took had none.
TEST_F(Simulation, KineticMcCountsTheJumpsWithoutABarrier) {
  const auto model_filename =
      pred::WriteModel(directory_ / "model_with_large_energy_changes.json",
                       pred::MakeSyntheticModelWithScaledEnergies(config_, pred::TestElementSet(), 600));
  std::string warnings;
  {
    const WorkingDirectory working_directory(directory_ / "kmc_without_barrier");
    const CapturedErrorOutput captured;
    // 300 steps, each of them logged, at a temperature at which the vacancy climbs about as often as it descends.
    KineticMcFirstSerial(config_, 1, 100, 300, 0, 0, 0, 0, 20000, pred::TestElementSet(), model_filename, "", false,
                         false, 41)
        .Simulate();
    warnings = captured.GetString();
  }
  const auto rows = ReadLog(directory_ / "kmc_without_barrier" / "kmc_log.txt");
  ASSERT_EQ(rows.size(), 301u);    // steps 0 to 300
  // columns: steps, time, temperature, energy, Ea, dE, ...
  size_t num_without_forward_barrier = 0;
  size_t num_without_backward_barrier = 0;
  std::string first_step;
  for (const auto &row: rows) {
    const double barrier = std::stod(row.at(4));
    const double energy_change = std::stod(row.at(5));
    ASSERT_GE(barrier, std::max(energy_change, 0.0)) << "step " << row.at(0);
    ASSERT_TRUE(std::isfinite(std::stod(row.at(1)))) << "step " << row.at(0);
    if (barrier == 0.0) {
      ++num_without_forward_barrier;
    } else if (barrier == energy_change) {
      ++num_without_backward_barrier;
    } else {
      continue;
    }
    if (first_step.empty()) {
      first_step = row.at(0);
    }
  }
  EXPECT_GE(num_without_forward_barrier, 10u);
  EXPECT_GE(num_without_backward_barrier, 10u);
  // The first one is reported when it happens, because a long run is stopped from outside and never gets to the
  // summary at its end.
  ASSERT_FALSE(first_step.empty());
  const std::string first_warning = "Warning: step " + first_step + ": so far 1 of " +
      std::to_string(std::stoull(first_step) + 1) + " jumps had no barrier in the quartic barrier model";
  EXPECT_EQ(warnings.rfind(first_warning, 0), 0u) << warnings;
  const std::string summary = "Warning: in this run " +
      std::to_string(num_without_forward_barrier + num_without_backward_barrier) +
      " of 301 jumps had no barrier in the quartic barrier model (" + std::to_string(num_without_forward_barrier) +
      " forward, " + std::to_string(num_without_backward_barrier) + " backward)";
  EXPECT_NE(warnings.find(summary), std::string::npos) << warnings;
}

// The warnings are the only output about jumps without a barrier, so a run inside the range of the quartic form
// writes what it wrote before they existed.
TEST_F(Simulation, KineticMcWritesNoWarningWhileEveryJumpHasABarrier) {
  std::string warnings;
  {
    const CapturedErrorOutput captured;
    RunKineticMc("kmc_with_barriers", 51);
    warnings = captured.GetString();
  }
  EXPECT_EQ(warnings, "");
  for (const auto &row: ReadLog(directory_ / "kmc_with_barriers" / "kmc_log.txt")) {
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

#include "VacancyMigrationPredictorQuartic.h"
#include "cfg/test/TemporaryPath.hpp"
#include "pred/test/SyntheticModel.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pred {
namespace {
// An FCC lattice of the given number of unit cells in a box with the given edges. The sites are numbered by position,
// as in a simulation.
cfg::Config MakeFccInBox(const Factor_t &factors, const Matrix_d &basis) {
  const auto fcc = cfg::GenerateFCC(factors, Element("Al"));
  cfg::Config config(basis, fcc.GetLatticeVector(), fcc.GetAtomVector(), fcc.GetMapShiftList(), true);
  config.ReassignLatticeVector();
  return config;
}
}    // namespace

// The barrier model needs a cubic box of at least four unit cells along every edge. Any other box used to be reported
// as a coefficient file of the wrong size, or as neighbourhoods that are not arranged alike.
class BarrierModelBox : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {
    temporary_directory_ = std::make_unique<test::TemporaryDirectory>();
    model_filename_ = WriteModel(temporary_directory_->GetPath() / "model.json",
                                 MakeSyntheticModel(cfg::GenerateFCC({4, 4, 4}, Element("Al")), TestElementSet()));
  }

  static void TearDownTestSuite() {
    temporary_directory_.reset();
  }

  // The error the predictor reports for this configuration, or nothing if it can be constructed.
  static std::string GetError(const cfg::Config &config, const std::string &model_filename = model_filename_) {
    try {
      const VacancyMigrationPredictorQuartic predictor(model_filename, config, TestElementSet());
    } catch (const std::runtime_error &error) {
      return error.what();
    }
    return "";
  }

  static void ExpectBoxError(const std::string &error, const std::string &box) {
    EXPECT_NE(error.find("The box of this configuration has " + box), std::string::npos) << error;
    EXPECT_NE(error.find("The barrier model needs a complete FCC lattice in a cubic box"), std::string::npos) << error;
  }

  static inline std::unique_ptr<test::TemporaryDirectory> temporary_directory_{};
  static inline std::string model_filename_{};
};

TEST_F(BarrierModelBox, CubesOfFourUnitCellsOrMoreAreAccepted) {
  for (size_t factor = 4; factor <= 6; ++factor) {
    EXPECT_EQ(GetError(cfg::GenerateFCC({factor, factor, factor}, Element("Al"))), "") << "factor " << factor;
  }
}

TEST_F(BarrierModelBox, CubeOfThreeUnitCellsIsRejected) {
  ExpectBoxError(GetError(cfg::GenerateFCC({3, 3, 3}, Element("Al"))),
                 "3 x 3 x 3 unit cells and the edges (12.138 0 0) (0 12.138 0) (0 0 12.138) Angstrom.");
}

// Which error such a box used to get depended on the first jump of the configuration, and with it on the numbers of the
// sites: a simulation numbers them by position, GenerateFCC() cell by cell.
TEST_F(BarrierModelBox, BoxWithUnequalEdgesIsRejected) {
  for (const Factor_t &factors: {Factor_t{6, 6, 7}, Factor_t{4, 4, 8}, Factor_t{8, 8, 12}, Factor_t{5, 6, 7},
                                 Factor_t{4, 8, 8}, Factor_t{3, 4, 4}}) {
    const auto box = std::to_string(factors[0]) + " x " + std::to_string(factors[1]) + " x " +
        std::to_string(factors[2]) + " unit cells";
    auto config = cfg::GenerateFCC(factors, Element("Al"));
    ExpectBoxError(GetError(config), box);
    config.ReassignLatticeVector();
    ExpectBoxError(GetError(config), box);
  }
}

// The other ways not to be a cube along x, y and z: an edge that is slightly too long, with which the encodings cannot
// even be built, a tilted edge, a cube that is turned, and one whose edges come in another order.
TEST_F(BarrierModelBox, BoxThatIsNotACubeAlongTheAxesIsRejected) {
  const double edge = 6 * constants::kLatticeConstant;
  const std::vector<std::pair<Matrix_d, std::string>> boxes{
      {{{{edge, 0, 0}, {0, edge, 0}, {0, 0, 1.001 * edge}}},
       "5.998 x 5.998 x 6.004 unit cells and the edges (24.276 0 0) (0 24.276 0) (0 0 24.3003) Angstrom."},
      {{{{edge, 0, 0}, {0.001 * edge, edge, 0}, {0, 0, edge}}},
       "6 x 6 x 6 unit cells and the edges (24.276 0 0) (0.024276 24.276 0) (0 0 24.276) Angstrom."},
      {{{{0.8 * edge, 0.6 * edge, 0}, {-0.6 * edge, 0.8 * edge, 0}, {0, 0, edge}}},
       "6 x 6 x 6 unit cells and the edges (19.4208 14.5656 0) (-14.5656 19.4208 0) (0 0 24.276) Angstrom."},
      {{{{0, edge, 0}, {0, 0, edge}, {edge, 0, 0}}},
       "6 x 6 x 6 unit cells and the edges (0 24.276 0) (0 0 24.276) (24.276 0 0) Angstrom."}};
  for (const auto &[basis, box]: boxes) {
    ExpectBoxError(GetError(MakeFccInBox({6, 6, 6}, basis)), box);
  }
}

TEST_F(BarrierModelBox, LatticeConstantAndRoundingOfACubeDoNotMatter) {
  for (const double lattice_constant: {4.0, 4.05, 4.1}) {
    const double edge = 4 * lattice_constant;
    EXPECT_EQ(GetError(MakeFccInBox({4, 4, 4}, {{{edge, 0, 0}, {0, edge, 0}, {0, 0, edge}}})), "")
        << "lattice constant " << lattice_constant;
  }
  const double edge = 4 * constants::kLatticeConstant;
  EXPECT_EQ(GetError(MakeFccInBox({4, 4, 4}, {{{edge, 1e-5, 0}, {0, edge + 1e-5, 0}, {0, 0, edge - 1e-5}}})), "");
}

TEST_F(BarrierModelBox, BoxIsCheckedBeforeTheCoefficientFileIsRead) {
  const auto missing_filename = (temporary_directory_->GetPath() / "missing.json").string();
  ExpectBoxError(GetError(cfg::GenerateFCC({6, 6, 7}, Element("Al")), missing_filename), "6 x 6 x 7 unit cells");
  EXPECT_EQ(GetError(cfg::GenerateFCC({6, 6, 6}, Element("Al")), missing_filename), "Cannot open " + missing_filename);
}

// Here the box is as it should be, and the comparison of the neighbourhoods of the jumps finds the mistake.
TEST_F(BarrierModelBox, LatticeWithAMissingSiteIsRejected) {
  const auto config = cfg::GenerateFCC({6, 6, 6}, Element("Al"));
  auto lattice_vector = config.GetLatticeVector();
  auto atom_vector = config.GetAtomVector();
  auto map_shift_list = config.GetMapShiftList();
  lattice_vector.pop_back();
  atom_vector.pop_back();
  map_shift_list.pop_back();
  const auto error = GetError({config.GetBasis(), lattice_vector, atom_vector, map_shift_list, true});
  EXPECT_NE(error.find("not arranged alike"), std::string::npos) << error;
  EXPECT_NE(error.find("The barrier model needs a complete FCC lattice in a cubic box"), std::string::npos) << error;
}
}    // namespace pred

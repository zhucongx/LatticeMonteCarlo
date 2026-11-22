#include "Config.h"

#include <gtest/gtest.h>

namespace cfg {
namespace {
Matrix_d IdentityBasis() {
  return {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
}

Config MakeLinearConfig(double lhs_x, double rhs_x) {
  auto basis = IdentityBasis();
  std::vector<Lattice> lattice_vector{Lattice(0, Vector_d{lhs_x, 0.0, 0.0}), Lattice(1, Vector_d{rhs_x, 0.0, 0.0})};
  std::vector<Atom> atom_vector{Atom(0, Element(ElementName::Al)), Atom(1, Element(ElementName::Mg))};
  std::vector<Vector_i> map_shift_list(2, Vector_i{0, 0, 0});
  return Config(basis, lattice_vector, atom_vector, map_shift_list, false);
}
}  // namespace

TEST(ConfigJump, HalfCellDisplacementDoesNotWrap) {
  auto config = MakeLinearConfig(0.0, 0.5);
  config.AtomJump({0, 1});

  const auto &map_shift_list = config.GetMapShiftList();
  EXPECT_EQ(map_shift_list[0][0], 0);
  EXPECT_EQ(map_shift_list[1][0], 0);
  EXPECT_DOUBLE_EQ(config.GetUnwrappedCartesianPositionOfAtom(0)[0], 0.5);
  EXPECT_DOUBLE_EQ(config.GetUnwrappedCartesianPositionOfAtom(1)[0], 0.0);
}

TEST(ConfigJump, WrapsAcrossPositiveBoundary) {
  auto config = MakeLinearConfig(0.99, 0.01);
  config.AtomJump({0, 1});

  const auto &map_shift_list = config.GetMapShiftList();
  EXPECT_EQ(map_shift_list[0][0], 1);
  EXPECT_EQ(map_shift_list[1][0], -1);
  // 0.99 - 1 is not the double closest to -0.01, so compare with a tolerance.
  EXPECT_NEAR(config.GetUnwrappedCartesianPositionOfAtom(0)[0], 1.01, 1e-12);
  EXPECT_NEAR(config.GetUnwrappedCartesianPositionOfAtom(1)[0], -0.01, 1e-12);
}

TEST(ConfigJump, WrapsFromNegativeStart) {
  auto config = MakeLinearConfig(-0.49, 0.49);
  config.LatticeJump({0, 1});

  // The constructor wraps the site at -0.49 to 0.51 and records that image in the map shift of its atom. The two
  // sites are then 0.02 apart, so the jump itself crosses no boundary and atom 1 keeps a map shift of zero.
  const auto &map_shift_list = config.GetMapShiftList();
  EXPECT_EQ(map_shift_list[0][0], -1);
  EXPECT_EQ(map_shift_list[1][0], 0);
  EXPECT_NEAR(config.GetUnwrappedCartesianPositionOfAtom(0)[0], -0.51, 1e-12);
  EXPECT_NEAR(config.GetUnwrappedCartesianPositionOfAtom(1)[0], 0.51, 1e-12);
}
}  // namespace cfg

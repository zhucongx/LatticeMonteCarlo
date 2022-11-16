#include <iostream>
#include <gtest/gtest.h>

#include "Constants.hpp"
#include "VectorMatrix.hpp"
#include "Lattice.hpp"
#include "Atom.hpp"
#include "LatticeCluster.hpp"
#include "Config.h"
#include "ClusterExpansionPeriodic.h"
#include "ClusterExpansionPredictorMm2.h"

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  cfg::Config a = cfg::Config::ReadCfg("test.cfg");
  auto lattice_id = a.GetVacancyLatticeIndex();
  EXPECT_EQ(lattice_id, 18);
  auto atom_id = a.GetLatticeToAtomHashmap().at(lattice_id);
  EXPECT_EQ(atom_id, 18);

  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
  EXPECT_EQ(7 * 6, 42);
  EXPECT_EQ(7 * 6, 42);
  EXPECT_EQ(7 * 6, 42);
}
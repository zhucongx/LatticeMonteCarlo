#include "Config.h"

#include <gtest/gtest.h>
#include <set>

namespace cfg {
// Supercells with fewer than three cutoff lengths per edge used to list the same neighbour several times.
TEST(ConfigNeighbors, EverySiteHasDistinctNeighborsInSmallSupercells) {
  for (size_t factor = 3; factor <= 5; ++factor) {
    const auto config = GenerateFCC({factor, factor, factor}, Element("Al"));
    for (size_t lattice_id = 0; lattice_id < config.GetNumAtoms(); ++lattice_id) {
      const auto &first = config.GetFirstNeighborsAdjacencyList()[lattice_id];
      const auto &second = config.GetSecondNeighborsAdjacencyList()[lattice_id];
      const auto &third = config.GetThirdNeighborsAdjacencyList()[lattice_id];
      ASSERT_EQ(first.size(), constants::kNumFirstNearestNeighbors) << "factor " << factor;
      ASSERT_EQ(second.size(), constants::kNumSecondNearestNeighbors) << "factor " << factor;
      ASSERT_EQ(third.size(), constants::kNumThirdNearestNeighbors) << "factor " << factor;
      ASSERT_EQ(std::set<size_t>(first.begin(), first.end()).size(), first.size()) << "factor " << factor;
      ASSERT_EQ(std::set<size_t>(second.begin(), second.end()).size(), second.size()) << "factor " << factor;
      ASSERT_EQ(std::set<size_t>(third.begin(), third.end()).size(), third.size()) << "factor " << factor;
    }
  }
}

TEST(ConfigNeighbors, NeighborListsAreSymmetric) {
  const auto config = GenerateFCC({3, 4, 5}, Element("Al"));
  for (size_t lattice_id = 0; lattice_id < config.GetNumAtoms(); ++lattice_id) {
    for (const auto neighbor_id: config.GetFirstNeighborsAdjacencyList()[lattice_id]) {
      EXPECT_EQ(config.FindDistanceLabelBetweenLattice(neighbor_id, lattice_id), 1);
    }
    for (const auto neighbor_id: config.GetThirdNeighborsAdjacencyList()[lattice_id]) {
      EXPECT_EQ(config.FindDistanceLabelBetweenLattice(neighbor_id, lattice_id), 3);
    }
  }
}

// Below two cutoff lengths a site sees two periodic images of the same neighbour.
TEST(ConfigNeighbors, BoxShorterThanTwoCutoffsIsRejected) {
  EXPECT_THROW(GenerateFCC({2, 2, 2}, Element("Al")), std::runtime_error);
  EXPECT_THROW(GenerateFCC({4, 4, 2}, Element("Al")), std::runtime_error);
}
}    // namespace cfg

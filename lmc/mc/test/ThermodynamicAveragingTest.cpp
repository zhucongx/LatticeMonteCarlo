#include "ThermodynamicAveraging.h"

#include <gtest/gtest.h>

namespace mc {
// The samples come from a trajectory that already visits states with their Boltzmann weight, so the average of
// the window is the plain mean. It used to be weighted by exp(-beta * E) a second time.
TEST(ThermodynamicAveraging, IsThePlainMeanOfTheWindow) {
  ThermodynamicAveraging averaging(3);
  averaging.AddEnergy(-1.0);
  EXPECT_DOUBLE_EQ(averaging.GetThermodynamicAverage(), -1.0);
  averaging.AddEnergy(-2.0);
  EXPECT_DOUBLE_EQ(averaging.GetThermodynamicAverage(), -1.5);
  averaging.AddEnergy(-6.0);
  EXPECT_DOUBLE_EQ(averaging.GetThermodynamicAverage(), -3.0);
}

TEST(ThermodynamicAveraging, DropsTheOldestSampleWhenTheWindowIsFull) {
  ThermodynamicAveraging averaging(3);
  for (const double energy: {-1.0, -2.0, -6.0, -10.0}) {
    averaging.AddEnergy(energy);
  }
  EXPECT_DOUBLE_EQ(averaging.GetThermodynamicAverage(), (-2.0 - 6.0 - 10.0) / 3.0);
}

TEST(ThermodynamicAveraging, WindowOfSizeZeroIsDisabled) {
  ThermodynamicAveraging averaging(0);
  averaging.AddEnergy(-5.0);
  EXPECT_DOUBLE_EQ(averaging.GetThermodynamicAverage(), 0.0);
}
}    // namespace mc

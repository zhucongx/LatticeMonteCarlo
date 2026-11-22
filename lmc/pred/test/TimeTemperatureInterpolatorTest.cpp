#include "TimeTemperatureInterpolator.h"
#include "cfg/test/TemporaryPath.hpp"

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

namespace pred {
namespace {
// Writes the text to a table file of its own and reads it.
TimeTemperatureInterpolator Read(const std::string &content) {
  const test::TemporaryFile file(content, ".dat");
  return TimeTemperatureInterpolator(file.GetName());
}
}    // namespace

TEST(TimeTemperatureInterpolator, InterpolatesLinearlyAndClampsOutsideTheTable) {
  const auto interpolator = Read("0 450\n10 350\n20 300\n");
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(0), 450);
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(5), 400);
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(15), 325);
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(-1), 450);
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(1000), 300);
}

// The last point used to be dropped when the file did not end with a newline.
TEST(TimeTemperatureInterpolator, KeepsTheLastPointOfAFileWithoutFinalNewline) {
  const auto interpolator = Read("0 450\n10 300");
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(10), 300);
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(5), 375);
}

// The reader used to search for a line that starts with the character '0' and never returned without one.
TEST(TimeTemperatureInterpolator, AcceptsHeaderCommentsBlanksAndAnyStartTime) {
  const auto interpolator = Read("time temperature\n# seconds Kelvin\n\n  1e-3\t450  # start\n2e-3 350\r\n");
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(1e-3), 450);
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(1.5e-3), 400);
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(1), 350);
}

TEST(TimeTemperatureInterpolator, SortsThePoints) {
  const auto interpolator = Read("10 300\n0 400\n");
  EXPECT_DOUBLE_EQ(interpolator.GetTemperature(5), 350);
}

TEST(TimeTemperatureInterpolator, RejectsUnusableTables) {
  EXPECT_THROW(Read(""), std::runtime_error);
  EXPECT_THROW(Read("time temperature\n"), std::runtime_error);
  EXPECT_THROW(Read("0 450\n10 -5\n"), std::runtime_error);
  EXPECT_THROW(Read("0 450\n0 300\n"), std::range_error);
  EXPECT_THROW(TimeTemperatureInterpolator(std::string("/no/such/directory/table.dat")), std::runtime_error);
}

TEST(TimeTemperatureInterpolator, EmptyFileNameMeansNoTable) {
  const std::string no_table;
  EXPECT_NO_THROW(TimeTemperatureInterpolator{no_table});
}
}    // namespace pred

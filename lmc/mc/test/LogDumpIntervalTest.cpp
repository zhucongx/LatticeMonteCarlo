#include "McAbstract.h"

#include <gtest/gtest.h>

namespace mc {
TEST(LogDumpInterval, GrowsByAFactorOfTenWithEveryDecadeOfSteps) {
  const unsigned long long int log_dump_steps = 100000;
  // Step 0 is where the previous floating point formula converted -1.0 to an unsigned integer.
  EXPECT_EQ(GetLogDumpInterval(0, log_dump_steps), 1u);
  EXPECT_EQ(GetLogDumpInterval(98, log_dump_steps), 1u);
  EXPECT_EQ(GetLogDumpInterval(99, log_dump_steps), 10u);
  EXPECT_EQ(GetLogDumpInterval(998, log_dump_steps), 10u);
  EXPECT_EQ(GetLogDumpInterval(999, log_dump_steps), 100u);
  EXPECT_EQ(GetLogDumpInterval(9999, log_dump_steps), 1000u);
  EXPECT_EQ(GetLogDumpInterval(99999, log_dump_steps), 10000u);
}

TEST(LogDumpInterval, NeverExceedsLogDumpSteps) {
  EXPECT_EQ(GetLogDumpInterval(999, 50), 50u);
  EXPECT_EQ(GetLogDumpInterval(2000, 250), 100u);
  EXPECT_EQ(GetLogDumpInterval(2501, 250), 250u);
  EXPECT_EQ(GetLogDumpInterval(1000000000000ULL, 100000), 100000u);
  EXPECT_EQ(GetLogDumpInterval(12345, 1), 1u);
}

TEST(LogDumpInterval, IsAPowerOfTenOrLogDumpSteps) {
  for (const unsigned long long int log_dump_steps: {1ULL, 7ULL, 100ULL, 250ULL, 100000ULL}) {
    for (unsigned long long int steps = 0; steps < 300000; ++steps) {
      const auto interval = GetLogDumpInterval(steps, log_dump_steps);
      ASSERT_GE(interval, 1u);
      ASSERT_LE(interval, log_dump_steps);
      if (interval != log_dump_steps) {
        auto power_of_ten = interval;
        while (power_of_ten % 10 == 0) {
          power_of_ten /= 10;
        }
        ASSERT_EQ(power_of_ten, 1u) << "steps " << steps << ", log_dump_steps " << log_dump_steps;
      }
    }
  }
}
}    // namespace mc

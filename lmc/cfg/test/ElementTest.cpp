#include "Element.hpp"

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

TEST(Element, KnownNamesRoundTrip) {
  for (const std::string name: {"X", "Al", "Mg", "Zn", "Cu", "Sn", "pAl", "pMg", "pZn", "pCu", "pSn"}) {
    EXPECT_EQ(Element(name).GetString(), name);
  }
}

// An unknown name used to become the vacancy "X", which turned a typo into an extra vacancy.
TEST(Element, UnknownNameIsRejected) {
  EXPECT_THROW(Element("Si"), std::invalid_argument);
  EXPECT_THROW(Element("al"), std::invalid_argument);
  EXPECT_THROW(Element(""), std::invalid_argument);
}

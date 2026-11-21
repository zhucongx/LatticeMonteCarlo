#include "EnergyUtility.h"
#include "pred/test/SyntheticModel.hpp"

#include <cstdint>
#include <gtest/gtest.h>
#include <random>
#include <stdexcept>
#include <vector>

namespace pred {
// OneHotEncoder replaces GetOneHotParametersFromMap() in the barrier model. The coefficients of the model were fitted
// to the encoding by element names, so the encoder has to reproduce it value for value.
TEST(OneHotEncoder, ReproducesTheEncodingByElementNames) {
  const auto config = MakeTestAlloy();
  const auto &element_set = TestElementSet();
  const std::vector<Element> elements(element_set.begin(), element_set.end());
  const auto encode_hashmap = GetOneHotEncodeHashmap(element_set);
  // The neighbourhood of a jump without the two sites of the jump.
  const size_t num_sites = constants::kNumThirdNearestSetSizeOfPair - 2;
  std::mt19937_64 generator(7);
  for (const auto &mapping: {GetAverageClusterParametersMappingMMM(config),
                             GetAverageClusterParametersMappingMM2(config)}) {
    const OneHotEncoder encoder(element_set, mapping);
    ASSERT_EQ(encoder.Size(), GetOneHotEncodeSize(element_set.size(), mapping));
    for (int trial = 0; trial < 50; ++trial) {
      std::vector<Element> site_elements;
      std::vector<uint8_t> site_element_indices;
      for (size_t i = 0; i < num_sites; ++i) {
        // Later trials are dominated by one element, as a dilute alloy is.
        const auto element = trial < 25 || generator() % 10 == 0 ? elements[generator() % elements.size()]
                                                                 : elements.front();
        site_elements.push_back(element);
        site_element_indices.push_back(encoder.GetElementIndex(element));
      }
      std::vector<double> expected;
      GetOneHotParametersFromMap(site_elements, encode_hashmap, element_set.size(), mapping, expected);
      std::vector<double> encode(encoder.Size(), -1.0);
      encoder.Encode(site_element_indices.data(), encode.data());
      ASSERT_EQ(encode, expected) << "trial " << trial;
    }
  }
}

TEST(OneHotEncoder, NumbersTheElementsInTheOrderOfTheElementSet) {
  const OneHotEncoder encoder(TestElementSet(), GetAverageClusterParametersMappingMMM(MakeTestAlloy()));
  EXPECT_EQ(encoder.GetElementIndex(Element("Al")), 0);
  EXPECT_EQ(encoder.GetElementIndex(Element("Mg")), 1);
  EXPECT_EQ(encoder.GetElementIndex(Element("Zn")), 2);
  // The encoding has no entry for the vacancy or for an element outside the element set.
  EXPECT_THROW((void) encoder.GetElementIndex(Element("X")), std::out_of_range);
  EXPECT_THROW((void) encoder.GetElementIndex(Element("Cu")), std::out_of_range);
}

TEST(OneHotEncoder, AveragesOverTheClustersOfAGroup) {
  // One group of three single sites, one group of two ordered pairs and one group of two pairs of equivalent sites.
  const OneHotEncoder encoder(TestElementSet(), {{{0}, {1}, {2}}, {{0, 1}, {1, 2}}, {{SIZE_MAX, 0, 1}, {SIZE_MAX, 2, 1}}});
  ASSERT_EQ(encoder.Size(), 3u + 9u + 6u);
  const std::vector<uint8_t> site_element_indices{0, 2, 0};    // Al Zn Al
  std::vector<double> encode(encoder.Size(), -1.0);
  encoder.Encode(site_element_indices.data(), encode.data());
  std::vector<double> expected(encoder.Size(), 0.0);
  expected[0] = 2.0 / 3.0;         // Al
  expected[2] = 1.0 / 3.0;         // Zn
  expected[3 + 0 * 3 + 2] = 0.5;   // Al, Zn in this order
  expected[3 + 2 * 3 + 0] = 0.5;   // Zn, Al in this order
  expected[12 + 2] = 1.0;          // Al and Zn in either order: the third of Al-Al, Al-Mg, Al-Zn, Mg-Mg, Mg-Zn, Zn-Zn
  EXPECT_EQ(encode, expected);
}

TEST(OneHotEncoder, RejectsMappingsItIsNotDefinedFor) {
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{0, 1, 2}}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{SIZE_MAX, 0, 1}, {0, 1}}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{0, 1}, {2}}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{300}}}), std::invalid_argument);
}
}    // namespace pred

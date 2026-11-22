#include "EnergyUtility.h"
#include "pred/test/SyntheticModel.hpp"

#include <algorithm>
#include <cmath>
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

// The barrier model does not build the encoding. It adds one weight per cluster, which has to be the dot product of
// the weights with the encoding.
TEST(OneHotEncoder, DotIsTheDotProductWithTheEncoding) {
  const auto config = MakeTestAlloy();
  const auto &element_set = TestElementSet();
  const size_t num_sites = constants::kNumThirdNearestSetSizeOfPair - 2;
  std::mt19937_64 generator(11);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  for (const auto &mapping: {GetAverageClusterParametersMappingMMM(config),
                             GetAverageClusterParametersMappingMM2(config)}) {
    const OneHotEncoder encoder(element_set, mapping);
    Eigen::VectorXd weights(static_cast<Eigen::Index>(encoder.Size()));
    for (Eigen::Index i = 0; i < weights.size(); ++i) {
      weights[i] = uniform(generator);
    }
    const auto scaled_weights = encoder.ScaleWeights(weights);
    for (int trial = 0; trial < 20; ++trial) {
      std::vector<uint8_t> site_element_indices(num_sites);
      for (auto &index: site_element_indices) {
        index = static_cast<uint8_t>(generator() % element_set.size());
      }
      std::vector<double> encode(encoder.Size());
      encoder.Encode(site_element_indices.data(), encode.data());
      const Eigen::Map<const Eigen::VectorXd> encode_vector(encode.data(), static_cast<Eigen::Index>(encode.size()));
      EXPECT_NEAR(encoder.Dot(site_element_indices.data(), scaled_weights.data()), weights.dot(encode_vector), 1e-12);
    }
    EXPECT_THROW((void) encoder.ScaleWeights(Eigen::VectorXd::Zero(3)), std::invalid_argument);
  }
}

// The frame of a jump is fixed up to a half turn about the jump direction. The mm2 encoding of the neighbourhood must
// not notice that half turn, otherwise the barrier of a jump would depend on the order of a neighbour list.
TEST(OneHotEncoder, Mm2EncodingDoesNotDependOnAHalfTurnAboutTheJumpDirection) {
  const auto config = MakeTestAlloy();
  const auto &element_set = TestElementSet();
  const std::pair<size_t, size_t> pair{0, config.GetFirstNeighborsAdjacencyList()[0][0]};
  const auto state_list = GetSortedLatticeVectorStateOfPair(config, pair);
  ASSERT_EQ(state_list.size(), constants::kNumThirdNearestSetSizeOfPair);
  // Entry k of the state list after half a turn about the x axis of the frame, which is the jump direction.
  std::vector<size_t> half_turn(state_list.size(), state_list.size());
  for (size_t k = 0; k < state_list.size(); ++k) {
    const auto &position = state_list[k].GetRelativePosition();
    for (size_t l = 0; l < state_list.size(); ++l) {
      const auto &image = state_list[l].GetRelativePosition();
      if (std::abs(image[0] - position[0]) < kEpsilon && std::abs(image[1] - (1 - position[1])) < kEpsilon &&
          std::abs(image[2] - (1 - position[2])) < kEpsilon) {
        half_turn[k] = l;
      }
    }
    ASSERT_LT(half_turn[k], state_list.size()) << "site " << k << " has no image";
  }
  const OneHotEncoder encoder(element_set, GetAverageClusterParametersMappingMM2(config));
  std::mt19937_64 generator(13);
  size_t num_moved_sites = 0;
  for (const auto &mm2_list: {GetSymmetricallySortedLatticeVectorMM2(config, pair),
                              GetSymmetricallySortedLatticeVectorMM2(config, {pair.second, pair.first})}) {
    for (int trial = 0; trial < 20; ++trial) {
      std::vector<uint8_t> state_element_indices(state_list.size());
      for (auto &index: state_element_indices) {
        index = static_cast<uint8_t>(generator() % element_set.size());
      }
      std::vector<uint8_t> elements;
      std::vector<uint8_t> turned_elements;
      for (const auto &lattice: mm2_list) {
        const auto it = std::find_if(state_list.begin(), state_list.end(), [&lattice](const cfg::Lattice &state_lattice) {
          return state_lattice.GetId() == lattice.GetId();
        });
        const auto position = static_cast<size_t>(it - state_list.begin());
        ASSERT_LT(position, state_list.size());
        if (half_turn[position] != position) {
          ++num_moved_sites;
        }
        elements.push_back(state_element_indices[position]);
        turned_elements.push_back(state_element_indices[half_turn[position]]);
      }
      std::vector<double> encode(encoder.Size());
      std::vector<double> turned_encode(encoder.Size());
      encoder.Encode(elements.data(), encode.data());
      encoder.Encode(turned_elements.data(), turned_encode.data());
      ASSERT_EQ(encode, turned_encode);
    }
  }
  EXPECT_GT(num_moved_sites, 0u);
}

TEST(OneHotEncoder, RejectsMappingsItIsNotDefinedFor) {
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{0, 1, 2}}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{SIZE_MAX, 0, 1}, {0, 1}}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{0, 1}, {2}}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{}}), std::invalid_argument);
  EXPECT_THROW(OneHotEncoder(TestElementSet(), {{{300}}}), std::invalid_argument);
}
}    // namespace pred

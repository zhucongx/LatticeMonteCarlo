#include "EnergyUtility.h"

#include <algorithm>
#include <gtest/gtest.h>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace pred {
namespace {
// The indexer the energy models build: every cluster of the elements and at most one vacancy.
ClusterIndexer MakeIndexer(std::set<Element> element_set) {
  element_set.emplace(ElementName::X);
  const auto cluster_hashmap = InitializeClusterHashMap(element_set);
  const std::map<cfg::ElementCluster, int> ordered(cluster_hashmap.begin(), cluster_hashmap.end());
  return ClusterIndexer(ordered, std::vector<double>(ordered.size(), 1.0));
}
}    // namespace

// The energy models count clusters through the lookup by element codes. It has to name the cluster that the lookup
// by cfg::ElementCluster names, whatever the order of the elements.
TEST(ClusterIndexer, LookupByCodesAgreesWithLookupByCluster) {
  const auto indexer = MakeIndexer({Element("Al"), Element("Mg"), Element("Zn")});
  ASSERT_GT(indexer.Size(), 0u);
  size_t num_orders = 0;
  for (const auto &cluster: indexer.GetClusters()) {
    const auto &elements = cluster.GetElementVector();
    std::vector<size_t> order(elements.size());
    std::iota(order.begin(), order.end(), 0);
    do {
      std::vector<uint8_t> codes;
      for (const auto position: order) {
        codes.push_back(indexer.GetElementCode(elements[position]));
      }
      EXPECT_EQ(indexer.GetIndex(static_cast<size_t>(cluster.GetLabel()), codes.data(), codes.size()),
                indexer.GetIndex(cluster));
      ++num_orders;
    } while (std::next_permutation(order.begin(), order.end()));
  }
  EXPECT_GT(num_orders, indexer.Size());
}

TEST(ClusterIndexer, DifferentElementsHaveDifferentCodes) {
  const auto indexer = MakeIndexer({Element("Al"), Element("Mg"), Element("Zn")});
  std::set<uint8_t> codes;
  for (const auto &name: {"Al", "Mg", "Zn", "X"}) {
    codes.insert(indexer.GetElementCode(Element(name)));
  }
  EXPECT_EQ(codes.size(), 4u);
  // Elements that the model does not know share one further code.
  EXPECT_EQ(indexer.GetElementCode(Element("Cu")), indexer.GetElementCode(Element("Sn")));
  EXPECT_EQ(codes.count(indexer.GetElementCode(Element("Cu"))), 0u);
}

TEST(ClusterIndexer, RejectsClustersThatAreNotPartOfTheModel) {
  const auto indexer = MakeIndexer({Element("Al"), Element("Mg")});
  const auto aluminium = indexer.GetElementCode(Element("Al"));
  const auto vacancy = indexer.GetElementCode(Element("X"));
  const auto foreign = indexer.GetElementCode(Element("Cu"));

  const uint8_t one_atom[] = {aluminium};
  const uint8_t valid_pair[] = {vacancy, aluminium};
  EXPECT_NO_THROW((void) indexer.GetIndex(0, one_atom, 1));
  EXPECT_NO_THROW((void) indexer.GetIndex(1, valid_pair, 2));

  // The model has no terms for two vacancies in one cluster or for an element outside the element set.
  const uint8_t two_vacancies[] = {vacancy, vacancy};
  const uint8_t foreign_atom[] = {foreign};
  const uint8_t foreign_pair[] = {aluminium, foreign};
  EXPECT_THROW((void) indexer.GetIndex(1, two_vacancies, 2), std::out_of_range);
  EXPECT_THROW((void) indexer.GetIndex(0, foreign_atom, 1), std::out_of_range);
  EXPECT_THROW((void) indexer.GetIndex(2, foreign_pair, 2), std::out_of_range);
  // A label of pairs asked for one element, and a label that does not exist.
  EXPECT_THROW((void) indexer.GetIndex(1, one_atom, 1), std::out_of_range);
  EXPECT_THROW((void) indexer.GetIndex(99, one_atom, 1), std::out_of_range);
  EXPECT_THROW((void) ClusterIndexer().GetIndex(0, one_atom, 1), std::out_of_range);
}

TEST(FlatClusterMapping, KeepsLabelsSizesAndSiteIndices) {
  const std::vector<std::vector<std::vector<size_t>>> mapping{{{5}, {7}}, {}, {{1, 2}, {3, 4}, {0, 59}}};
  const auto flat = FlattenClusterMapping(mapping);
  ASSERT_EQ(flat.groups.size(), 2u);
  EXPECT_EQ(flat.groups[0].label, 0u);
  EXPECT_EQ(flat.groups[0].cluster_size, 1u);
  EXPECT_EQ(flat.groups[0].num_clusters, 2u);
  EXPECT_EQ(flat.groups[1].label, 2u);
  EXPECT_EQ(flat.groups[1].cluster_size, 2u);
  EXPECT_EQ(flat.groups[1].num_clusters, 3u);
  EXPECT_EQ(flat.groups[1].offset, 2u);
  EXPECT_EQ(flat.site_indices, (std::vector<uint8_t>{5, 7, 1, 2, 3, 4, 0, 59}));

  EXPECT_THROW(FlattenClusterMapping({{{1, 2}, {3}}}), std::invalid_argument);
  EXPECT_THROW(FlattenClusterMapping({{{1, 2, 3, 4}}}), std::invalid_argument);
  EXPECT_THROW(FlattenClusterMapping({{{300}}}), std::invalid_argument);
}
}    // namespace pred

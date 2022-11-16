#include "ShortRangeOrder.h"

#include <utility>

namespace ansys {
ShortRangeOrder::ShortRangeOrder(cfg::Config config, std::set<Element> element_set)
    : config_(std::move(config)), element_set_(std::move(element_set)) {}
std::map<std::string, double> ShortRangeOrder::FindWarrenCowley() const {
  std::map<std::string, double> res;
  for (const auto &element1: element_set_) {
    for (const auto &element2: element_set_) {
      res[element1.GetString() + "-" + element2.GetString()] = 0;
    }
  }
  const auto &element_list_map = config_.GetElementAtomIdVectorMap();
  std::map<Element, double> concentration;
  std::map<Element, double> count;
  for (const auto type: element_set_) {
    count[type] = static_cast<double>(element_list_map.at(type).size());
    concentration[type] = static_cast<double>(element_list_map.at(type).size())
        / static_cast<double>(config_.GetNumAtoms());
  }

  for (const auto &atom1: config_.GetAtomVector()) {
    const auto element1 = atom1.GetElement();
    if (element1 == ElementName::X) {
      continue;
    }
    const auto &neighbor_list1 = config_.GetFirstNeighborsAtomIdVectorOfAtom(atom1.GetId());
    for (const auto &element2: element_set_) {
      double bond_count = 0;
      for (auto neighbor_index: neighbor_list1) {
        if (config_.GetAtomVector()[neighbor_index].GetElement() == element2) {
          bond_count += 1;
        }
      }
      double pij = bond_count / 12;
      double aij = (pij - concentration[element2])
          / (static_cast<double>(atom1.GetElement() == element2) - concentration[element2]);

      res[element1.GetString() + "-" + element2.GetString()] += aij / count[element1];
    }
  }
  // for (const auto &type1: element_set_) {
  //   for (const auto &type2: element_set_) {
  //     std::cout << std::string(type1).append("-").append(type2) << ' '
  //               << res[std::string(type1).append("-").append(type2)] << '\n';
  //   }
  // }

  return res;
}

} // ansys
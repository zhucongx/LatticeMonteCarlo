#include "ShortRangeOrder.h"

#include <utility>

namespace ansys {
ShortRangeOrder::ShortRangeOrder(cfg::Config config, std::set<Element> element_set)
    : config_(std::move(config)), element_set_(std::move(element_set)) {}
std::map<std::string, double> ShortRangeOrder::FindWarrenCowley(const size_t shell_number) const {
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
  size_t num_bonds;
  switch (shell_number) {
    case 1:num_bonds = constants::kNumFirstNearestNeighbors;
      break;
    case 2:num_bonds = constants::kNumSecondNearestNeighbors;
      break;
    case 3:num_bonds = constants::kNumThirdNearestNeighbors;
      break;
    default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
  }
  for (const auto &element1: element_set_) {
    size_t num_all_bonds = element_list_map.at(element1).size() * num_bonds;
    std::map<Element, size_t> ct_this_pair_map{};
    for (const auto &atom_id1: element_list_map.at(element1)) {
      std::vector<size_t> neighbor_list;
      switch (shell_number) {
        case 1:neighbor_list = config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id1);
          break;
        case 2:neighbor_list = config_.GetSecondNeighborsAtomIdVectorOfAtom(atom_id1);
          break;
        case 3:neighbor_list = config_.GetThirdNeighborsAtomIdVectorOfAtom(atom_id1);
          break;
      }
      for (const auto &atom_id2: neighbor_list) {
        ct_this_pair_map[config_.GetElementAtAtomId(atom_id2)]++;
      }
    }
    for (auto [element2, ct_this_pair]: ct_this_pair_map) {
      std::string key = element1.GetString() + "-" + element2.GetString();
      double pij = static_cast<double>(ct_this_pair) / static_cast<double>(num_all_bonds);
      res[key] = (pij - concentration[element2])
          / ((static_cast<double>(element1 == element2) - concentration[element2]));
    }
  }
  return res;
}

} // ansys
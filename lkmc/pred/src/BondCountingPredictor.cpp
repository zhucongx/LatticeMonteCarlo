#include "BondCountingPredictor.h"
namespace pred {
std::unordered_map<ElementBond, size_t, boost::hash<ElementBond> > InitializeHashMap(
    const std::set<Element> &type_set) {
  std::unordered_map<ElementBond, size_t, boost::hash<ElementBond> > bond_count_hashmap;
  for (size_t label = 1; label <= 7; ++label) {
    for (const auto &element1: type_set) {
      for (const auto &element2: type_set) {
        bond_count_hashmap[ElementBond(label, element1, element2)] = 0;
      }
    }
  }
  return bond_count_hashmap;
}
std::vector<double> GetBondChange(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    std::unordered_map<ElementBond, size_t, boost::hash<ElementBond> > initialized_hashmap) {
  auto vacancy_lattice_id = lattice_id_jump_pair.first;
  auto migration_atom_jump_id = lattice_id_jump_pair.second;
  const Element type1 = config.GetElementAtLatticeId(migration_atom_jump_id);
  // plus new bonds
  for (auto lattice2_id: config.GetFirstNeighborsAdjacencyList()[vacancy_lattice_id]) {
    if (lattice2_id == migration_atom_jump_id) { continue; }
    initialized_hashmap[ElementBond{1, type1, config.GetElementAtLatticeId(lattice2_id)}]++;
  }
  for (auto lattice2_id: config.GetSecondNeighborsAdjacencyList()[vacancy_lattice_id]) {
    initialized_hashmap[ElementBond{2, type1, config.GetElementAtLatticeId(lattice2_id)}]++;
  }
  for (auto lattice2_id: config.GetThirdNeighborsAdjacencyList()[vacancy_lattice_id]) {
    initialized_hashmap[ElementBond{3, type1, config.GetElementAtLatticeId(lattice2_id)}]++;
  }
  for (auto lattice2_id: config.GetFourthNeighborsAdjacencyList()[vacancy_lattice_id]) {
    initialized_hashmap[ElementBond{4, type1, config.GetElementAtLatticeId(lattice2_id)}]++;
  }
  for (auto lattice2_id: config.GetFifthNeighborsAdjacencyList()[vacancy_lattice_id]) {
    initialized_hashmap[ElementBond{5, type1, config.GetElementAtLatticeId(lattice2_id)}]++;
  }
  for (auto lattice2_id: config.GetSixthNeighborsAdjacencyList()[vacancy_lattice_id]) {
    initialized_hashmap[ElementBond{6, type1, config.GetElementAtLatticeId(lattice2_id)}]++;
  }
  for (auto lattice2_id: config.GetSeventhNeighborsAdjacencyList()[vacancy_lattice_id]) {
    initialized_hashmap[ElementBond{7, type1, config.GetElementAtLatticeId(lattice2_id)}]++;
  }
  // minus old bonds
  for (auto lattice2_id: config.GetFirstNeighborsAdjacencyList()[migration_atom_jump_id]) {
    if (lattice2_id == vacancy_lattice_id) { continue; }
    initialized_hashmap[ElementBond{1, type1, config.GetElementAtLatticeId(lattice2_id)}]--;
  }
  for (auto lattice2_id: config.GetSecondNeighborsAdjacencyList()[migration_atom_jump_id]) {
    initialized_hashmap[ElementBond{2, type1, config.GetElementAtLatticeId(lattice2_id)}]--;
  }
  for (auto lattice2_id: config.GetThirdNeighborsAdjacencyList()[migration_atom_jump_id]) {
    initialized_hashmap[ElementBond{3, type1, config.GetElementAtLatticeId(lattice2_id)}]--;
  }
  for (auto lattice2_id: config.GetFourthNeighborsAdjacencyList()[migration_atom_jump_id]) {
    initialized_hashmap[ElementBond{4, type1, config.GetElementAtLatticeId(lattice2_id)}]--;
  }
  for (auto lattice2_id: config.GetFifthNeighborsAdjacencyList()[migration_atom_jump_id]) {
    initialized_hashmap[ElementBond{5, type1, config.GetElementAtLatticeId(lattice2_id)}]--;
  }
  for (auto lattice2_id: config.GetSixthNeighborsAdjacencyList()[migration_atom_jump_id]) {
    initialized_hashmap[ElementBond{6, type1, config.GetElementAtLatticeId(lattice2_id)}]--;
  }
  for (auto lattice2_id: config.GetSeventhNeighborsAdjacencyList()[migration_atom_jump_id]) {
    initialized_hashmap[ElementBond{7, type1, config.GetElementAtLatticeId(lattice2_id)}]--;
  }

  std::map<ElementBond, int> ordered(initialized_hashmap.begin(), initialized_hashmap.end());
  std::vector<double> res;
  res.reserve(ordered.size());
  for (const auto &bond_count: ordered) {
    res.push_back(bond_count.second);
  }
  return res;
}
} // namespace pred
#ifndef LKMC_LKMC_PRED_INCLUDE_BONDCOUNTINGPREDICTOR_H_
#define LKMC_LKMC_PRED_INCLUDE_BONDCOUNTINGPREDICTOR_H_
#include "ElementBond.hpp"
#include "Config.h"
namespace pred {
std::unordered_map<ElementBond, size_t, boost::hash<ElementBond> > InitializeHashMap(
    const std::set<Element> &type_set);
std::vector<double> GetBondChange(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    std::unordered_map<ElementBond, size_t, boost::hash<ElementBond> > initialized_hashmap);
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_BONDCOUNTINGPREDICTOR_H_

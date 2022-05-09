
#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORUTILITY_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORUTILITY_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
namespace pred {
using Singlet_MMM_t = cfg::LatticeClusterMMM<1>;
using Pair_MMM_t = cfg::LatticeClusterMMM<2>;
using Singlet_MM2_t = cfg::LatticeClusterMM2<1>;
using Pair_MM2_t = cfg::LatticeClusterMM2<2>;

using Singlet_State_t = cfg::LatticeCluster<1>;
using Pair_State_t = cfg::LatticeCluster<2>;
using Triplet_State_t = cfg::LatticeCluster<3>;

std::unordered_map<std::string, std::vector<double> > GetOneHotEncodeHashmap(
    const std::set<Element> &type_set);
// Returns forward and backward sorted lattice lists
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMMM(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMM2(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMM2(
    const cfg::Config &config);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMMM(
    const cfg::Config &config);

std::unordered_map<
    cfg::ElementCluster, int, boost::hash<cfg::ElementCluster> > InitializeClusterHashMap(
    const std::set<Element> &type_set);
// Returns forward and backward sorted lattice lists
std::vector<cfg::Lattice> GetSortedLatticeVectorState(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t> > > GetClusterParametersMappingState(
    const cfg::Config &config);
std::vector<double> GetOneHotParametersFromMap(
    const std::vector<Element> &encode,
    const std::unordered_map<std::string, std::vector<double> > &one_hot_encode_hashmap,
    size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t> > > &cluster_mapping);
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORUTILITY_H_

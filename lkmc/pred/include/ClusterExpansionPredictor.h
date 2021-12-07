#ifndef LKMC_LKMC_PRED_INCLUDE_CLUSTEREXPANSIONPREDICTOR_H_
#define LKMC_LKMC_PRED_INCLUDE_CLUSTEREXPANSIONPREDICTOR_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
namespace pred {
std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMMM(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMMM(
    const cfg::Config &config);

// std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMM2(
//     const cfg::Config &config);
// std::vector<cfg::Lattice> GetSymmetricallySortedLatticeVectorMM2(
//     const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);

std::vector<cfg::Lattice> GetSortedLatticeVectorPeriodic(
    const cfg::Config &config, size_t lattice_id);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingPeriodic(
    const cfg::Config &config);

struct Parameters_E0 {
  std::vector<double> mu_x{};
  std::vector<double> sigma_x{};
  double mu_y{};
  double sigma_y{};
  std::vector<std::vector<double> > U{};
  std::vector<double> theta{};
};
struct Parameters_dE {
  std::vector<double> mu_x{};
  std::vector<double> sigma_x{};
  double mu_y{};
  double sigma_y{};
  std::vector<double> theta{};
};
class ClusterExpansionPredictor {
  public:
    ClusterExpansionPredictor(const std::string &predictor_filename,
                              const cfg::Config &reference_config,
                              const std::set<Element> &type_set);
    virtual ~ClusterExpansionPredictor();
    [[nodiscard]]  std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
    [[nodiscard]]  std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  protected:
    [[nodiscard]]  double GetE0(const cfg::Config &config,
                                const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                Element migration_element) const;
    [[nodiscard]]  double GetE(const cfg::Config &config,
                               size_t lattice_id,
                               Element migration_element) const;
  private:
    size_t type_size_;

    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;

    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_{};
    const std::vector<std::vector<std::vector<size_t> > > mapping_periodic_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;
    std::unordered_map<size_t,
                       std::vector<size_t> > site_cluster_periodic_hashmap_;

    std::unordered_map<Element,
                       Parameters_E0,
                       boost::hash<Element>> element_parameters_hashmap_;
    Parameters_dE parameters_dE_;
};

} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_CLUSTEREXPANSIONPREDICTOR_H_

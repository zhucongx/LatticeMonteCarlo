#ifndef LKMC_LKMC_PRED_INCLUDE_CLUSTEREXPANSIONPREDICTOR_H_
#define LKMC_LKMC_PRED_INCLUDE_CLUSTEREXPANSIONPREDICTOR_H_
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
namespace pred {
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMMM(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingMM2(
    const cfg::Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::vector<std::vector<std::vector<size_t> > > GetAverageClusterParametersMappingPeriodic(
    const cfg::Config &config, size_t lattice_id_jump_pair);
struct Element_Parameters {
  std::vector<double> mu_x{};
  std::vector<std::vector<double> > transform_matrix{};
  std::vector<double> theta{};
  double mean_y{};
};

class ClusterExpansionPredictor {
  public:
    ClusterExpansionPredictor(const std::string &predictor_filename,
                              const cfg::Config &reference_config,
                              const std::set<std::string> &type_set);
    virtual ~ClusterExpansionPredictor();

  private:
    std::unordered_map<std::string, Element_Parameters> element_parameters_hashmap_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > >
        site_bond_cluster_mmm_hashmap_{};
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > >
        site_bond_cluster_mm2_hashmap_{};
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > >
        site_bond_cluster_periodic_hashmap_{};
};

} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_CLUSTEREXPANSIONPREDICTOR_H_

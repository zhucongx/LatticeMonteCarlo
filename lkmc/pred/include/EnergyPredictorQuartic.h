#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "ElementCluster.hpp"
#include "EnergyPredictorUtility.h"

namespace pred {
struct ParametersQuartic {
  std::vector<double> mu_x_mmm{};
  std::vector<double> mu_x_mm2{};
  std::vector<double> sigma_x_mmm{};
  std::vector<double> sigma_x_mm2{};
  std::vector<std::vector<double> > U_mmm{};
  std::vector<std::vector<double> > U_mm2{};
  std::vector<double> theta_D{};
  std::vector<double> theta_Ks{};
  double mu_y_D{};
  double mu_y_Ks{};
  double sigma_y_D{};
  double sigma_y_Ks{};
};

class EnergyPredictorQuartic {
  public:
    EnergyPredictorQuartic(const std::string &predictor_filename,
                           const cfg::Config &reference_config,
                           std::set<Element> type_set);
    virtual ~EnergyPredictorQuartic();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
  protected:
    [[nodiscard]] double GetDe(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetKs(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] double GetD(const cfg::Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;


    const std::set<Element> type_set_;
    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_mm2_;
    const std::vector<std::vector<std::vector<size_t> > > mapping_state_;

    std::vector<double> base_theta_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mm2_hashmap_;
    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_state_hashmap_;

    std::unordered_map<Element,
                       ParametersQuartic,
                       boost::hash<Element>> element_parameters_hashmap_;

    std::unordered_map<cfg::ElementCluster, int,
                       boost::hash<cfg::ElementCluster> > initialized_cluster_hashmap_;

};

} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_

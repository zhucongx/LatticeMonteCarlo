#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "LatticeCluster.hpp"
namespace pred {
struct Parameters {
  std::vector<double> e0_mu_x{};
  std::vector<double> e0_sigma_x{};
  double e0_mu_y{};
  double e0_sigma_y{};
  std::vector<std::vector<double> > e0_U{};
  std::vector<double> e0_theta{};

  std::vector<double> dE_mu_x{};
  std::vector<double> dE_sigma_x{};
  double dE_mu_y{};
  double dE_sigma_y{};
  std::vector<std::vector<double> > dE_U{};
  std::vector<double> dE_theta{};
};

// class EnergyPredictorQuartic {
//     EnergyPredictorQuartic(const std::string &predictor_filename,
//                            const cfg::Config &reference_config,
//                            const std::set<Element> &type_set);
//     virtual ~EnergyPredictorQuartic();
//     [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
//         const cfg::Config &config,
//         const std::pair<size_t, size_t> &atom_id_jump_pair) const;
//   protected:
//
//     [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
//         const cfg::Config &config,
//         const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
//
//     [[nodiscard]] double GetA(const cfg::Config &config,
//                               const std::pair<size_t, size_t> &lattice_id_jump_pair,
//                               Element migration_element) const;
//     [[nodiscard]] double GetB(const cfg::Config &config,
//                               const std::pair<size_t, size_t> &lattice_id_jump_pair,
//                               Element migration_element) const;
//     [[nodiscard]] double GetC(const cfg::Config &config,
//                               const std::pair<size_t, size_t> &lattice_id_jump_pair,
//                               Element migration_element) const;
//
//   private:
//     const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
//     const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_;
//     const std::vector<std::vector<std::vector<size_t> > > mapping_mm2_;
//
//     std::unordered_map<std::pair<size_t, size_t>,
//                        std::vector<size_t>,
//                        boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;
//     std::unordered_map<std::pair<size_t, size_t>,
//                        std::vector<size_t>,
//                        boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mm2_hashmap_;
//     std::unordered_map<Element,
//                        Parameters,
//                        boost::hash<Element>> e0_element_parameters_hashmap_;
//
// };

} // namespace pred
#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_

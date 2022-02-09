#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#include <string>
#include <set>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "ClusterExpansionPredictor.h"
namespace pred {
struct ParametersE0 {
  std::vector<double> mu_x{};
  std::vector<double> sigma_x{};
  double mu_y{};
  double sigma_y{};
  std::vector<std::vector<double> > U{};
  std::vector<double> theta{};
};
struct ParametersDE {
  std::vector<double> mu_x{};
  std::vector<double> sigma_x{};
  double mu_y{};
  double sigma_y{};
  std::vector<std::vector<double> > U{};
  std::vector<double> theta{};
};
struct ParametersQuartic {
  std::vector<double> mu_x_mmm{};
  std::vector<double> mu_x_mm2{};
  std::vector<double> sigma_x_mmm{};
  std::vector<double> sigma_x_mm2{};
  std::vector<std::vector<double> > U_mmm{};
  std::vector<std::vector<double> > U_mm2{};
  std::vector<double> theta_a{};
  std::vector<double> theta_b{};
  std::vector<double> theta_c{};
  double mu_a{};
  double mu_b{};
  double mu_c{};
  double sigma_a{};
  double sigma_b{};
  double sigma_c{};
};
class EnergyPredictor {
  public:
    explicit EnergyPredictor(const std::set<Element> &type_set);
    virtual ~EnergyPredictor();
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromAtomIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const;
  protected:
    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const = 0;
    std::set<Element> type_set_;
    const std::unordered_map<std::string, std::vector<double> > one_hot_encode_hash_map_;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_
#include "EnergyPredictor.h"
namespace preds {
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
class EnergyPredictorQuartic : public EnergyPredictor{
  public:
    EnergyPredictorQuartic(const std::string &predictor_filename,
                                     const cfg::Config &reference_config,
                                     const std::set<Element> &type_set);
    ~EnergyPredictorQuartic() override;
  protected:
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const override;
    [[nodiscard]] std::pair<double,double> GetAAndC(const cfg::Config &config,
                                const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                Element migration_element) const;
    [[nodiscard]] double GetB(const cfg::Config &config,
                               const std::pair<size_t, size_t> &lattice_id_jump_pair,
                               Element migration_element) const;
  private:
    const std::vector<std::vector<std::vector<size_t> > > mapping_mmm_{};
    const std::vector<std::vector<std::vector<size_t> > > mapping_mm2_{};

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mmm_hashmap_;

    std::unordered_map<std::pair<size_t, size_t>,
                       std::vector<size_t>,
                       boost::hash<std::pair<size_t, size_t> > > site_bond_cluster_mm2_hashmap_;
    std::unordered_map<Element,
                       ParametersQuartic,
                       boost::hash<Element>> element_parameters_hashmap_;
};
} // namespace preds

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTIC_H_

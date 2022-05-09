#ifndef LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTICLRU_H_
#define LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTICLRU_H_
#include "EnergyPredictorQuartic.h"
namespace pred {

class EnergyPredictorQuarticLru : public EnergyPredictorQuartic {
  public:
    EnergyPredictorQuarticLru(const std::string &predictor_filename,
                              const cfg::Config &reference_config,
                              const std::set<Element> &type_set,
                              size_t cache_size);
    ~EnergyPredictorQuarticLru() override;
  private:
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const override;
    [[nodiscard]] size_t GetHashFromConfigAndLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

    void Add(size_t key, std::pair<double, double> value) const;
    mutable std::list<std::pair<size_t, std::pair<double, double> > > cache_list_{};
    mutable std::unordered_map<
        size_t,
        std::list<std::pair<size_t, std::pair<double, double> > >::iterator> hashmap_{};
    const size_t cache_size_;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_ENERGYPREDICTORQUARTICLRU_H_

#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
#include "LruCache.hpp"
#include "VacancyMigrationPredictorQuartic.h"

namespace pred {

class VacancyMigrationPredictorQuarticLru : public VacancyMigrationPredictorQuartic {
 public:
  VacancyMigrationPredictorQuarticLru(const std::string &predictor_filename,
                                      const cfg::Config &reference_config,
                                      const std::set<Element> &element_set,
                                      size_t cache_size);
  ~VacancyMigrationPredictorQuarticLru() override;
  [[nodiscard]] std::pair<double, double>
  GetBarrierAndDiffFromLatticeIdPair(const cfg::Config &config,
                                     const std::pair<size_t, size_t> &lattice_id_jump_pair) const override;

 private:
  [[nodiscard]] size_t GetHashFromConfigAndLatticeIdPair(const cfg::Config &config,
                                                         const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  mutable LruCache<size_t, std::pair<double, double>> lru_cache_;
};
}    // namespace pred

#endif    //LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_

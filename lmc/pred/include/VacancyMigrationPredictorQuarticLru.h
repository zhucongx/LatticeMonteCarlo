#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
#include "LruCache.hpp"
#include "VacancyMigrationPredictorQuartic.h"

#include <array>
#include <cstdint>
#include <cstring>

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
  // The elements on the sites around a jump, one byte per site in the order of the state list of the pair. The barrier
  // and the energy change are a function of exactly these elements, and the key holds them in full: two neighbourhoods
  // share a cache entry only if they are the same. (A hash of the elements, used as the key before, would hand out
  // the barrier of another neighbourhood whenever two hashes collide.)
  static constexpr size_t kMaxKeySites = 64;
  using CacheKey = std::array<uint8_t, kMaxKeySites>;
  struct CacheKeyHash {
    size_t operator()(const CacheKey &key) const noexcept {
      uint64_t hash = 0x9E3779B97F4A7C15ULL;
      for (size_t offset = 0; offset < key.size(); offset += sizeof(uint64_t)) {
        uint64_t word;
        std::memcpy(&word, key.data() + offset, sizeof(uint64_t));
        hash = (hash ^ word) * 0xFF51AFD7ED558CCDULL;
        hash ^= hash >> 32;
      }
      return static_cast<size_t>(hash);
    }
  };
  [[nodiscard]] CacheKey GetKeyFromConfigAndLatticeIdPair(const cfg::Config &config,
                                                          const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  mutable LruCache<CacheKey, std::pair<double, double>, CacheKeyHash> lru_cache_;
};
}    // namespace pred

#endif    //LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
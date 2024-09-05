#include "VacancyMigrationPredictorE0Lru.h"
namespace pred {

VacancyMigrationPredictorE0Lru::VacancyMigrationPredictorE0Lru(const std::string &predictor_filename,
                                                               const cfg::Config &reference_config,
                                                               const std::set<Element> &element_set,
                                                               size_t cache_size)
    : VacancyMigrationPredictorE0(predictor_filename, reference_config, element_set),
      lru_cache_(cache_size) {}
VacancyMigrationPredictorE0Lru::~VacancyMigrationPredictorE0Lru() = default;
std::pair<double, double> VacancyMigrationPredictorE0Lru::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto key = GetHashFromConfigAndLatticeIdPair(config, lattice_id_jump_pair);
  std::pair<double, double> value;
  if (lru_cache_.Get(key, value)) {
    return value;
  } else {
    value = VacancyMigrationPredictorE0::GetBarrierAndDiffFromLatticeIdPair(
        config, lattice_id_jump_pair);
    lru_cache_.Add(key, value);
  }
  return value;
}
size_t VacancyMigrationPredictorE0Lru::GetHashFromConfigAndLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto &lattice_id_list_state = site_bond_cluster_state_hashmap_.at(lattice_id_jump_pair);
  const auto &lattice_id_list_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);

  size_t seed = 0;
  for (size_t i = 0; i < constants::kNumThirdNearestSetSizeOfPair; i++) {
    boost::hash_combine(seed, config.GetAtomIdFromLatticeId(lattice_id_list_state[i]));
  }
  for (size_t i = 0; i < (constants::kNumThirdNearestSetSizeOfPair - 2); i++) {
    boost::hash_combine(seed, config.GetAtomIdFromLatticeId(lattice_id_list_mmm[i]));
  }

  return seed;
}
} // pred

#include "VacancyMigrationPredictorQuarticLru.h"
namespace pred {

VacancyMigrationPredictorQuarticLru::VacancyMigrationPredictorQuarticLru(const std::string &predictor_filename,
                                                                         const cfg::Config &reference_config,
                                                                         const std::set<Element> &element_set,
                                                                         const size_t cache_size)
    : VacancyMigrationPredictorQuartic(predictor_filename, reference_config, element_set),
      lru_cache_(cache_size) {}
VacancyMigrationPredictorQuarticLru::~VacancyMigrationPredictorQuarticLru() = default;
std::pair<double, double> VacancyMigrationPredictorQuarticLru::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto key = GetHashFromConfigAndLatticeIdPair(config, lattice_id_jump_pair);
  std::pair<double, double> value;
  if (lru_cache_.Get(key, value)) {
    return value;
  } else {
    value = VacancyMigrationPredictorQuartic::GetBarrierAndDiffFromLatticeIdPair(
        config, lattice_id_jump_pair);
    lru_cache_.Add(key, value);
  }
  return value;
}
size_t VacancyMigrationPredictorQuarticLru::GetHashFromConfigAndLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  // TODO(arch/perf): Flatten pair->neighbor tables into contiguous arrays (site*12+nn_idx)
  // and precompute element-type signatures so this hot path avoids unordered_map<pair>
  // lookups, .at() bounds checks, and repeated hash_combine.
  const auto flat_idx = GetPairFlatIndex(lattice_id_jump_pair);
  const auto &lattice_id_list_state = site_bond_cluster_state_flat_.at(flat_idx);
  const auto &lattice_id_list_mmm = site_bond_cluster_mmm_flat_.at(flat_idx);
  const auto &lattice_id_list_mm2 = site_bond_cluster_mm2_flat_.at(flat_idx);

  size_t seed = 0;
  // Use element-type signature instead of atom ids to define the key.
  // This increases cache hit rate and reduces key-construction overhead.
  for (size_t i = 0; i < constants::kNumThirdNearestSetSizeOfPair; i++) {
    boost::hash_combine(seed, hash_value(config.GetElementAtLatticeId(lattice_id_list_state[i])));
  }
  for (size_t i = 0; i < (constants::kNumThirdNearestSetSizeOfPair - 2); i++) {
    boost::hash_combine(seed, hash_value(config.GetElementAtLatticeId(lattice_id_list_mmm[i])));
  }
  for (size_t i = 0; i < (constants::kNumThirdNearestSetSizeOfPair - 2); i++) {
    boost::hash_combine(seed, hash_value(config.GetElementAtLatticeId(lattice_id_list_mm2[i])));
  }
  return seed;
}
} // pred

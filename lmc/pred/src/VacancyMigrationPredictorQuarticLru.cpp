#include "VacancyMigrationPredictorQuarticLru.h"
namespace pred {

VacancyMigrationPredictorQuarticLru::VacancyMigrationPredictorQuarticLru(const std::string &predictor_filename,
                                                                         const cfg::Config &reference_config,
                                                                         const std::set<Element> &element_set,
                                                                         size_t cache_size)
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
  // TODO(arch): Consider flattening the pair->neighbor lookup tables into
  // contiguous arrays (site*12 + nn_idx) to avoid unordered_map<pair> lookups
  // and bounds-checked .at() on hot paths.
  const auto &lattice_id_list_state = site_bond_cluster_state_hashmap_.at(lattice_id_jump_pair);
  const auto &lattice_id_list_mmm = site_bond_cluster_mmm_hashmap_.at(lattice_id_jump_pair);
  const auto &lattice_id_list_mm2 = site_bond_cluster_mm2_hashmap_.at(lattice_id_jump_pair);

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
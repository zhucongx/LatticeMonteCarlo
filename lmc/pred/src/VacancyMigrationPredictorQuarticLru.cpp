#include "VacancyMigrationPredictorQuarticLru.h"

#include <stdexcept>
#include <string>
namespace pred {

VacancyMigrationPredictorQuarticLru::VacancyMigrationPredictorQuarticLru(const std::string &predictor_filename,
                                                                         const cfg::Config &reference_config,
                                                                         const std::set<Element> &element_set,
                                                                         const size_t cache_size)
    : VacancyMigrationPredictorQuartic(predictor_filename, reference_config, element_set),
      lru_cache_(cache_size) {
  for (const auto &lattice_id_list_state: site_bond_cluster_state_flat_) {
    if (lattice_id_list_state.size() > kMaxKeySites) {
      throw std::runtime_error("The neighbourhood of a jump has " + std::to_string(lattice_id_list_state.size()) +
                               " sites, the barrier cache holds at most " + std::to_string(kMaxKeySites));
    }
  }
}
VacancyMigrationPredictorQuarticLru::~VacancyMigrationPredictorQuarticLru() = default;
std::pair<double, double> VacancyMigrationPredictorQuarticLru::GetBarrierAndDiffFromLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  const auto key = GetKeyFromConfigAndLatticeIdPair(config, lattice_id_jump_pair);
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
VacancyMigrationPredictorQuarticLru::CacheKey VacancyMigrationPredictorQuarticLru::GetKeyFromConfigAndLatticeIdPair(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  // The lists of the mmm and mm2 encodings hold the sites of the state list in another order, which is the same for
  // every pair (the constructor of the base class checks it). The elements in the order of the state list therefore
  // determine every input of the model.
  const auto &lattice_id_list_state = site_bond_cluster_state_flat_[GetPairFlatIndex(lattice_id_jump_pair)];
  CacheKey key{};
  for (size_t i = 0; i < lattice_id_list_state.size(); ++i) {
    key[i] = static_cast<uint8_t>(hash_value(config.GetElementAtLatticeId(lattice_id_list_state[i])));
  }
  return key;
}
} // pred

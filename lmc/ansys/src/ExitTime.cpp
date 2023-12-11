#include "ExitTime.h"
namespace ansys {

ExitTime::ExitTime(const cfg::Config &config,
                   const Element &solvent_element,
                   const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor,
                   const double temperature)
    : config_(config),
      solvent_element_(solvent_element),
      vacancy_migration_predictor_(vacancy_migration_predictor),
      beta_(1.0 / constants::kBoltzmann / temperature)
{
}

std::pair<std::vector<std::vector<double>>, std::vector<double>> ExitTime::GetBarrierListAndExitTime() const
{
  std::vector<std::vector<double>> barrier_lists{};
  auto this_config = config_;
  this_config.SetAtomElementTypeAtAtom(this_config.GetVacancyAtomId(), solvent_element_);

  for (size_t atom_id = 0; atom_id < config_.GetNumAtoms(); ++atom_id) {
    std::cout << "atom_id: " << atom_id << std::endl;
    const size_t lattice_id = this_config.GetLatticeIdFromAtomId(atom_id);
    const Element this_element = this_config.GetElementAtAtomId(atom_id);
    this_config.SetAtomElementTypeAtAtom(atom_id, Element(ElementName::X));

    std::vector<double> barrier_list{};
    barrier_list.reserve(constants::kNumFirstNearestNeighbors);
    for (auto neighbor_lattice_id: this_config.GetFirstNeighborsAdjacencyList()[lattice_id]) {
      barrier_list.push_back(vacancy_migration_predictor_
                                 .GetBarrierAndDiffFromLatticeIdPair(this_config, {lattice_id, neighbor_lattice_id})
                                 .first);
    }
    this_config.SetAtomElementTypeAtAtom(atom_id, this_element);

    std::sort(barrier_list.begin(), barrier_list.end());
    barrier_lists.push_back(barrier_list);
  }
  std::vector<double> exit_times{};
  exit_times.reserve(config_.GetNumAtoms());
  for (const auto &barriers: barrier_lists) {
    double total_rate_k = 0.0;
    for (const auto &barrier: barriers) {
      total_rate_k += std::exp(-barrier * beta_);
    }
    exit_times.push_back(1.0 / total_rate_k / constants::kPrefactor);
  }
  return std::make_pair(barrier_lists, exit_times);
}

}    // namespace ansys

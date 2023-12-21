#include "ExitTime.h"

namespace ansys {
ExitTime::ExitTime(const cfg::Config &config,
                   const Element &solvent_element,
                   const double temperature,
                   const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor,
                   const pred::EnergyChangePredictorPairSite &energy_change_predictor_site,
                   const std::map<Element, double> &chemical_potential)
    : config_(config),
      solvent_element_(solvent_element),
      beta_(1.0 / constants::kBoltzmann / temperature),
      vacancy_migration_predictor_(vacancy_migration_predictor),
      energy_change_predictor_pair_site_(energy_change_predictor_site),
      chemical_potential_(chemical_potential) {}

std::pair<std::vector<std::vector<double>>, std::vector<double>> ExitTime::GetBarrierListAndExitTime() const {
  std::vector<std::vector<double>> barrier_lists{};
  auto this_config = config_;
  this_config.SetAtomElementTypeAtAtom(this_config.GetVacancyAtomId(), solvent_element_);
  for (size_t atom_id = 0; atom_id < this_config.GetNumAtoms(); ++atom_id) {
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
  exit_times.reserve(this_config.GetNumAtoms());
  for (const auto &barriers: barrier_lists) {
    double total_rate_k = 0.0;
    for (const auto &barrier: barriers) {
      total_rate_k += std::exp(-barrier * beta_);
    }
    exit_times.push_back(1.0 / total_rate_k / constants::kPrefactor);
  }
  return std::make_pair(barrier_lists, exit_times);
}

std::map<Element, std::vector<double>> ExitTime::GetBindingEnergy() const {
  std::map<Element, std::vector<double>> binding_energies{};
  auto this_config = config_;
  this_config.SetAtomElementTypeAtAtom(this_config.GetVacancyAtomId(), solvent_element_);
  for (size_t atom_id = 0; atom_id < this_config.GetNumAtoms(); ++atom_id) {
    const Element this_element = this_config.GetElementAtAtomId(atom_id);
    for (const auto &[element, mu]: chemical_potential_) {
      if (element == Element("X")) {
        continue;
      }
      this_config.SetAtomElementTypeAtAtom(atom_id, element);
      const auto potential_change =
          chemical_potential_.at(Element("X")) - chemical_potential_.at(this_config.GetElementAtAtomId(atom_id));
      const auto energy_change =
          energy_change_predictor_pair_site_.GetDeFromAtomIdSite(this_config, atom_id, Element("X"));

      if (binding_energies.find(element) == binding_energies.end()) {
        binding_energies.insert({element, {}});
      }

      binding_energies.at(element).push_back(energy_change - potential_change);
    }
    this_config.SetAtomElementTypeAtAtom(atom_id, this_element);
  }
  return binding_energies;
}

// double ExitTime::GetLocalBindingEnergy() const {
//   const auto potential_change = chemical_potential_.at(Element("X"));
//
//   const auto energy_change =
//       energy_change_predictor_pair_site_.GetDeFromAtomIdSite(config_, config_.GetVacancyAtomId(), solvent_element_);
//
//   return -energy_change - potential_change;
// }

// std::vector<double> ExitTime::GetBindingEnergy() const {
//   std::vector<double> binding_energies{};
//   auto this_config = config_;
//   this_config.SetAtomElementTypeAtAtom(this_config.GetVacancyAtomId(), solvent_element_);
//   for (size_t atom_id = 0; atom_id < this_config.GetNumAtoms(); ++atom_id) {
//     const auto potential_change =
//         chemical_potential_.at(Element("X")) - chemical_potential_.at(this_config.GetElementAtAtomId(atom_id));
//
//     const auto energy_change =
//         energy_change_predictor_pair_site_.GetDeFromAtomIdSite(this_config, atom_id, Element("X"));
//
//     binding_energies.push_back(energy_change - potential_change);
//   }
//   return binding_energies;
// }
// std::vector<double> ExitTime::GetProfileEnergy() const {
//   std::vector<double> profile_energies{};
//   const auto vacancy_atom_id = config_.GetVacancyAtomId();
//
//   const auto potential_change = chemical_potential_.at(Element("X"));
//   const auto energy_change =
//       energy_change_predictor_pair_site_.GetDeFromAtomIdSite(config_, vacancy_atom_id, solvent_element_);
//   const auto energy_reference = -energy_change - potential_change;
//   for (size_t atom_id = 0; atom_id < config_.GetNumAtoms(); ++atom_id) {
//     if (atom_id == vacancy_atom_id) {
//       profile_energies.push_back(energy_reference);
//       continue;
//     }
//     profile_energies.push_back(
//         energy_reference + energy_change_predictor_pair_site_.GetDeFromAtomIdPair(config_, {atom_id, vacancy_atom_id}));
//   }
//   return profile_energies;
// }
}    // namespace ansys

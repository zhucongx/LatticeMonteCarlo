#include "ExitTime.h"

#include "Eigen/Dense"

namespace ansys {
ExitTime::ExitTime(const cfg::Config &config,
                   const Element &solvent_element,
                   std::set<Element> element_set,
                   const double temperature,
                   const pred::EnergyPredictor &energy_predictor_,
                   const pred::VacancyMigrationPredictorQuartic &vacancy_migration_predictor,
                   const pred::EnergyChangePredictorPairSite &energy_change_predictor_site,
                   const std::map<Element, double> &chemical_potential)
    : config_(config),
      solvent_element_(solvent_element),
      element_set_(std::move(element_set)),
      beta_(1.0 / constants::kBoltzmann / temperature),
      energy_predictor_(energy_predictor_),
      vacancy_migration_predictor_(vacancy_migration_predictor),
      energy_change_predictor_pair_site_(energy_change_predictor_site),
      chemical_potential_(chemical_potential) {}

void ExitTime::GetExitTimeInfo(nlohmann::json &frame_info,
                               std::map<std::string, cfg::Config::VectorVariant> &auxiliary_lists,
                               std::map<std::string, cfg::Config::ValueVariant> &global_list) const {
  const auto profile_energy_list = GetProfileEnergy();
  auxiliary_lists["profile_energy"] = profile_energy_list;

  const auto binding_energy = GetBindingEnergy();
  std::vector<double> binding_energy_list;
  binding_energy_list.reserve(config_.GetNumAtoms());
  for (size_t atom_id = 0; atom_id < config_.GetNumAtoms(); ++atom_id) {
    std::vector<double> element_energy_list;
    element_energy_list.reserve(element_set_.size());
    for (const auto &element: element_set_) {
      element_energy_list.push_back(binding_energy.at(element)[atom_id]);
    }
    binding_energy_list.push_back(*std::max_element(element_energy_list.begin(), element_energy_list.end()));
  }
  auxiliary_lists["binding_energy"] = binding_energy_list;

  std::map<std::string, double> vac_element_energy_list;
  for (const auto &element: element_set_) {
    const auto element_string = element.GetString();
    const auto vac_element_energy = binding_energy.at(element)[config_.GetVacancyAtomId()];
    vac_element_energy_list.emplace(element_string, vac_element_energy);
    auxiliary_lists["binding_energy_" + element_string] = binding_energy.at(element);
  }
  frame_info["vacancy_local_binding_energy"] = vac_element_energy_list;


  std::unordered_set<size_t> all_atom_id_set{};
  for (const auto &cluster_info: frame_info["clusters"]) {
    const auto atom_id_set = cluster_info["cluster_atom_id_list"].get<std::unordered_set<size_t>>();
    all_atom_id_set.insert(atom_id_set.begin(), atom_id_set.end());
  }
  const auto [pair_energy_map, neighbor_atom_id_lists, migration_barrier_lists, driving_force_lists] =
      GetJumpEnergetics(all_atom_id_set);
  // auxiliary_lists["neighbor_atom_id_lists"] = neighbor_atom_id_lists;
  // auxiliary_lists["migration_barrier_lists"] = migration_barrier_lists;
  // auxiliary_lists["driving_force_lists"] = driving_force_lists;

  for (auto &cluster_info: frame_info["clusters"]) {
    std::vector<double> cluster_binding_energy_list;
    std::vector<double> cluster_profile_energy_list;
    std::map<Element, std::vector<double>> cluster_binding_energy_list_map;
    for (const auto &element: element_set_) {
      if (element == Element("X")) {
        continue;
      }
      cluster_binding_energy_list_map.insert({element, {}});
    }
    for (const auto &atom_id: cluster_info["cluster_atom_id_list"]) {
      std::vector<double> element_energy_list;
      element_energy_list.reserve(element_set_.size());
      for (const auto &element: element_set_) {
        cluster_binding_energy_list_map.at(element).push_back(binding_energy.at(element)[atom_id]);
        element_energy_list.push_back(binding_energy.at(element)[atom_id]);
      }
      cluster_binding_energy_list.push_back(*std::max_element(element_energy_list.begin(), element_energy_list.end()));
      cluster_profile_energy_list.push_back(profile_energy_list[atom_id]);
    }
    cluster_info["vacancy_binding_energy"] =
        *std::min_element(cluster_binding_energy_list.begin(), cluster_binding_energy_list.end());
    cluster_info["vacancy_profile_energy"] =
        *std::min_element(cluster_profile_energy_list.begin(), cluster_profile_energy_list.end());

    const auto atom_id_set = cluster_info["cluster_atom_id_list"].get<std::unordered_set<size_t>>();
    std::unordered_set<size_t> atom_id_set_plus_nn{};
    for (const auto &atom_id: atom_id_set) {
      atom_id_set_plus_nn.insert(atom_id);
      for (const auto neighbor_atom_id: neighbor_atom_id_lists[atom_id]) {
        atom_id_set_plus_nn.insert(neighbor_atom_id);
      }
    }
    // for (const auto &atom_id: atom_id_set) {
    //   atom_id_set_plus_nn.insert(atom_id);
    //   size_t lattice_id = config_.GetLatticeIdFromAtomId(atom_id);
    //   for (const auto &neighbor_lattice_id: config_.GetFirstNeighborsAdjacencyList()[lattice_id]) {
    //     atom_id_set_plus_nn.insert(config_.GetAtomIdFromLatticeId(neighbor_lattice_id));
    //   }
    //   for (const auto &neighbor_lattice_id: config_.GetSecondNeighborsAdjacencyList()[lattice_id]) {
    //     atom_id_set_plus_nn.insert(config_.GetAtomIdFromLatticeId(neighbor_lattice_id));
    //   }
    //   for (const auto &neighbor_lattice_id: config_.GetThirdNeighborsAdjacencyList()[lattice_id]) {
    //     atom_id_set_plus_nn .insert(config_.GetAtomIdFromLatticeId(neighbor_lattice_id));
    //   }
    // }

    cluster_info["markov_escape_time"] =
        BuildMarkovChain(atom_id_set_plus_nn, neighbor_atom_id_lists, migration_barrier_lists, profile_energy_list);

    cluster_info["barriers"] = GetAverageBarriers(atom_id_set, pair_energy_map);
    for (const auto &element: element_set_) {
      if (element == Element("X")) {
        continue;
      }
      cluster_info["vacancy_binding_energy_" + element.GetString()] = *std::min_element(
          cluster_binding_energy_list_map.at(element).begin(), cluster_binding_energy_list_map.at(element).end());
    }
    // cluster_info.erase("cluster_atom_id_list");
  }

  // auto [barrier_lists, average_barriers, exit_times] = GetBarrierListAndExitTime();
  // auxiliary_lists["barrier_lists"] = barrier_lists;
  // auxiliary_lists["average_barriers"] = average_barriers;
  // auxiliary_lists["exit_times"] = exit_times;
}

double ExitTime::BuildMarkovChain(const std::unordered_set<size_t> &atom_id_set,
                                  const std::vector<std::vector<size_t>> &neighbor_atom_id_lists,
                                  const std::vector<std::vector<double>> &migration_barrier_lists,
                                  const std::vector<double> &base_energy_list) const {

  const int transient_size = static_cast<int>(atom_id_set.size());

  Eigen::MatrixXd transient_matrix = Eigen::MatrixXd::Zero(transient_size, transient_size);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(transient_size, transient_size);
  Eigen::VectorXd probability_vector = Eigen::VectorXd::Zero(transient_size);
  Eigen::VectorXd tau_vector = Eigen::VectorXd::Zero(transient_size);

  std::unordered_map<size_t, int> atom_index_map;
  for (int index = 0; index < transient_size; ++index) {
    auto atom_id = *std::next(atom_id_set.begin(), index);
    atom_index_map.insert({atom_id,index});
  }

  for (const auto &atom_id: atom_id_set) {
    int index = atom_index_map[atom_id];
    probability_vector(index) = std::exp(-base_energy_list.at(atom_id) * beta_);

    std::vector<double> rates;
    double total_rate = 0.0;

    for (size_t j = 0; j < constants::kNumFirstNearestNeighbors; ++j) {
      double barrier = migration_barrier_lists.at(atom_id).at(j);
      double rate = constants::kPrefactor * std::exp(-barrier * beta_);
      rates.push_back(rate);
      total_rate += rate;
    }
    tau_vector(index) = 1.0 / total_rate;

    // Populate the Markov matrix
    for (size_t j = 0; j < constants::kNumFirstNearestNeighbors; ++j) {
      size_t neighbor_atom_id = neighbor_atom_id_lists.at(atom_id).at(j);
      auto it = atom_index_map.find(neighbor_atom_id);
      if (it == atom_index_map.end()) {
        // Jump to absorbing state
        continue;
      }
      auto neighbor_idx = it->second;
      transient_matrix(index, neighbor_idx) = rates.at(j) / total_rate;
    }
  }
  probability_vector /= probability_vector.sum();
  const double escape_time = probability_vector.transpose() * (I - transient_matrix).inverse() * tau_vector;
  return escape_time;
}

std::tuple<std::vector<std::vector<double>>, std::vector<double>, std::vector<double>>
ExitTime::GetBarrierListAndExitTime() const {
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
  std::vector<double> average_barriers{};
  std::vector<double> exit_times{};
  exit_times.reserve(this_config.GetNumAtoms());
  for (const auto &barriers: barrier_lists) {
    double sum_barrier = 0.0;
    double total_rate_k = 0.0;
    for (const auto &barrier: barriers) {
      sum_barrier += barrier;
      total_rate_k += std::exp(-barrier * beta_);
    }
    average_barriers.push_back(sum_barrier / static_cast<double>(barriers.size()));
    exit_times.push_back(1.0 / total_rate_k / constants::kPrefactor);
  }
  return std::make_tuple(barrier_lists, average_barriers, exit_times);
}

std::tuple<
    std::unordered_map<std::pair<size_t, size_t>, std::pair<double, double>, boost::hash<std::pair<size_t, size_t>>>,
    std::vector<std::vector<size_t>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<double>>>
ExitTime::GetJumpEnergetics(const std::unordered_set<size_t> &atom_id_set) const {
  std::vector<std::vector<size_t>> neighbor_atom_id_lists(
      config_.GetNumAtoms(), std::vector<size_t>(constants::kNumFirstNearestNeighbors, config_.GetNumAtoms()));
  std::vector<std::vector<double>> migration_barrier_lists(
      config_.GetNumAtoms(), std::vector<double>(constants::kNumFirstNearestNeighbors, nan("")));
  std::vector<std::vector<double>> driving_force_lists(
      config_.GetNumAtoms(), std::vector<double>(constants::kNumFirstNearestNeighbors, nan("")));

  std::unordered_set<size_t> lattice_id_set{};
  std::unordered_set<size_t> lattice_id_set_plus_nn{};
  for (const auto &atom_id: atom_id_set) {
    const auto lattice_id = config_.GetLatticeIdFromAtomId(atom_id);
    lattice_id_set.insert(lattice_id);
    lattice_id_set_plus_nn.insert(lattice_id);
    for (const auto &neighbor_lattice_id: config_.GetFirstNeighborsAdjacencyList()[lattice_id]) {
      lattice_id_set_plus_nn.insert(neighbor_lattice_id);
    }
    // for (const auto &neighbor_lattice_id: config_.GetSecondNeighborsAdjacencyList()[lattice_id]) {
    //   lattice_id_set_plus_nn.insert(neighbor_lattice_id);
    // }
    // for (const auto &neighbor_lattice_id: config_.GetThirdNeighborsAdjacencyList()[lattice_id]) {
    //   lattice_id_set_plus_nn.insert(neighbor_lattice_id);
    // }
  }

  std::unordered_map<std::pair<size_t, size_t>, std::pair<double, double>, boost::hash<std::pair<size_t, size_t>>>
      energy_barrier_map{};
  auto this_config = config_;
  this_config.SetAtomElementTypeAtAtom(this_config.GetVacancyAtomId(), solvent_element_);

  for (const auto lattice_id: lattice_id_set_plus_nn) {
    const auto atom_id = this_config.GetAtomIdFromLatticeId(lattice_id);
    std::vector<size_t> neighbor_atom_id_list{};
    std::vector<double> migration_barrier_list{};
    std::vector<double> driving_force_list{};

    const Element this_element = this_config.GetElementAtLatticeId(lattice_id);
    this_config.SetAtomElementTypeAtLattice(lattice_id, Element(ElementName::X));
    for (const auto neighbor_lattice_id: this_config.GetFirstNeighborsAdjacencyList().at(lattice_id)) {
      neighbor_atom_id_list.emplace_back(this_config.GetAtomIdFromLatticeId(neighbor_lattice_id));
      const auto jump_pair = std::make_pair(lattice_id, neighbor_lattice_id);
      const auto energetic_pair =
          vacancy_migration_predictor_.GetBarrierAndDiffFromLatticeIdPair(this_config, jump_pair);
      energy_barrier_map.insert({jump_pair, energetic_pair});
      migration_barrier_list.push_back(energetic_pair.first);
      driving_force_list.push_back(energetic_pair.second);
    }
    neighbor_atom_id_lists[atom_id] = neighbor_atom_id_list;
    migration_barrier_lists[atom_id] = migration_barrier_list;
    driving_force_lists[atom_id] = driving_force_list;
    this_config.SetAtomElementTypeAtLattice(lattice_id, this_element);
  }
  return std::make_tuple(energy_barrier_map, neighbor_atom_id_lists, migration_barrier_lists, driving_force_lists);
}

std::vector<double>
ExitTime::GetAverageBarriers(const std::unordered_set<size_t> &atom_id_set,
                             const std::unordered_map<std::pair<size_t, size_t>,
                                                      std::pair<double, double>,
                                                      boost::hash<std::pair<size_t, size_t>>> &pair_energy_map) const {
  std::unordered_set<size_t> lattice_id_set{};
  std::unordered_set<size_t> lattice_id_set_plus_nn{};
  for (const auto &atom_id: atom_id_set) {
    const auto lattice_id = config_.GetLatticeIdFromAtomId(atom_id);
    lattice_id_set.insert(lattice_id);
    lattice_id_set_plus_nn.insert(lattice_id);
    for (const auto &neighbor_lattice_id: config_.GetFirstNeighborsAdjacencyList()[lattice_id]) {
      lattice_id_set_plus_nn.insert(neighbor_lattice_id);
    }
  }
  std::vector<double> barrier_list_in{};
  std::vector<double> barrier_list_to{};
  std::vector<double> barrier_list_on{};
  std::vector<double> barrier_list_off{};
  auto this_config = config_;
  this_config.SetAtomElementTypeAtAtom(this_config.GetVacancyAtomId(), solvent_element_);

  for (const auto lattice_id: lattice_id_set_plus_nn) {
    const bool lattice_in_set = lattice_id_set.find(lattice_id) != lattice_id_set.end();

    const Element this_element = this_config.GetElementAtLatticeId(lattice_id);
    this_config.SetAtomElementTypeAtLattice(lattice_id, Element(ElementName::X));

    for (const auto neighbor_lattice_id: this_config.GetFirstNeighborsAdjacencyList().at(lattice_id)) {
      const bool neighbor_in_set = lattice_id_set.find(neighbor_lattice_id) != lattice_id_set.end();
      const bool neighbor_in_nn = lattice_id_set_plus_nn.find(neighbor_lattice_id) != lattice_id_set_plus_nn.end();

      std::vector<double> *barrier_list = nullptr;
      if (lattice_in_set && neighbor_in_set) {
        barrier_list = &barrier_list_in;
      } else if (lattice_in_set && !neighbor_in_set && neighbor_in_nn) {
        barrier_list = &barrier_list_to;
      } else if (!lattice_in_set && !neighbor_in_set && neighbor_in_nn) {
        barrier_list = &barrier_list_on;
      } else if (!lattice_in_set && !neighbor_in_nn) {
        barrier_list = &barrier_list_off;
      }
      if (barrier_list) {
        const auto jump_pair = std::make_pair(lattice_id, neighbor_lattice_id);
        double Ea, dE;
        if (pair_energy_map.find(jump_pair) == pair_energy_map.end()) {
          std::tie(Ea, dE) = vacancy_migration_predictor_.GetBarrierAndDiffFromLatticeIdPair(
              this_config, {lattice_id, neighbor_lattice_id});
        } else {
          std::tie(Ea, dE) = pair_energy_map.at({lattice_id, neighbor_lattice_id});
        }
        barrier_list->push_back(Ea - dE / 2);
      }
    }
    this_config.SetAtomElementTypeAtLattice(lattice_id, this_element);
  }

  auto compute_average = [](const std::vector<double> &barrier_list) -> double {
    if (barrier_list.empty()) {
      return nan("");
    }
    return std::accumulate(barrier_list.begin(), barrier_list.end(), 0.0) / static_cast<double>(barrier_list.size());
  };

  auto compute_std = [](const std::vector<double> &barrier_list, const double mean) -> double {
    if (barrier_list.empty()) {
      return nan("");
    }
    double sum = 0;
    for (const auto &barrier: barrier_list) {
      sum += (barrier - mean) * (barrier - mean);
    }
    return std::sqrt(sum / static_cast<double>(barrier_list.size()));
  };

  const double mean_in = compute_average(barrier_list_in);
  const double mean_to = compute_average(barrier_list_to);
  const double mean_on = compute_average(barrier_list_on);
  const double mean_off = compute_average(barrier_list_off);

  const double std_in = compute_std(barrier_list_in, mean_in);
  const double std_to = compute_std(barrier_list_to, mean_to);
  const double std_on = compute_std(barrier_list_on, mean_on);
  const double std_off = compute_std(barrier_list_off, mean_off);


  return {mean_in, mean_to, mean_on, mean_off, std_in, std_to, std_on, std_off};
}

std::map<Element, std::vector<double>> ExitTime::GetBindingEnergy() const {
  std::map<Element, std::vector<double>> binding_energies{};
  for (const auto &[element, mu]: chemical_potential_) {
    if (element == Element("X")) {
      continue;
    }
    binding_energies.insert({element, {}});
  }

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
// std::vector<double> ExitTime::GetSwapEnergy() const {
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

std::vector<double> ExitTime::GetProfileEnergy() const {
  std::vector<double> formation_energies{};

  const auto vacancy_atom_id = config_.GetVacancyAtomId();
  auto this_config = config_;
  this_config.SetAtomElementTypeAtAtom(vacancy_atom_id, solvent_element_);

  // const double  total_energy = energy_predictor_.GetEnergy(this_config);
  // const auto real_chemical_potential = energy_predictor_.GetRealChemicalPotential(solvent_element_);

  // double formation_energy_reference = total_energy;
  // for (const auto &[element, mu]: real_chemical_potential) {
  //   if (element == Element("X")) {
  //     continue;
  //   }
  //   formation_energy_reference -= mu * static_cast<double>(element_count_map.at(element));
  // }
  // std::unordered_set<size_t> neighbor_atom_set = this_config.GetNeighborsAtomIdSetOfAtom(vacancy_atom_id);

  for (size_t atom_id = 0; atom_id < this_config.GetNumAtoms(); ++atom_id) {
    if (atom_id == vacancy_atom_id) {
      std::vector<double> try_list;
      for (const auto &element: element_set_) {
        if (element == Element("X")) {
          continue;
        }
        this_config.SetAtomElementTypeAtAtom(vacancy_atom_id, element);
        const double energy_change = energy_change_predictor_pair_site_.GetDeFromAtomIdSite(this_config, atom_id, Element("X"));
        try_list.push_back(energy_change + chemical_potential_.at(element) - chemical_potential_.at(Element("X")));
      }
      formation_energies.push_back(*std::max_element(try_list.begin(), try_list.end()));
      this_config.SetAtomElementTypeAtAtom(vacancy_atom_id, solvent_element_);
      continue;
    }
    const Element this_element = this_config.GetElementAtAtomId(atom_id);
    const double energy_change = energy_change_predictor_pair_site_.GetDeFromAtomIdSite(this_config, atom_id, Element("X"));
    formation_energies.push_back(energy_change + chemical_potential_.at(this_element) - chemical_potential_.at(Element("X")));
  }
  return formation_energies;
}
}    // namespace ansys

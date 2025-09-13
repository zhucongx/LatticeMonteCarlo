#include "SimulatedAnnealing.h"

#include "EnergyPredictor.h"

#include <chrono>
#include <utility>

namespace mc {

static std::set<Element> GetElementSetFromSolventAndSolute(Element solvent_element,
                                                           const std::map<Element, size_t> &solute_atom_count) {
  std::set<Element> element_set;
  element_set.insert(solvent_element);
  for (const auto &solute_element_count: solute_atom_count) {
    element_set.insert(solute_element_count.first);
  }
  return element_set;
}

static std::vector<size_t> GetSoluteAtomIdVector(const cfg::Config &config,
                                                 const std::map<Element, size_t> &solute_atom_count) {
  std::vector<size_t> solute_atom_id_vector;
  for (const auto &atom: config.GetAtomVector()) {
    if (solute_atom_count.find(atom.GetElement()) != solute_atom_count.end()) {
      solute_atom_id_vector.push_back(atom.GetId());
    }
  }
  return solute_atom_id_vector;
}

SimulatedAnnealing::SimulatedAnnealing(const Factor_t &factors,
                                       Element solvent_element,
                                       const std::map<Element, size_t> &solute_atom_count,
                                       const unsigned long long int log_dump_steps,
                                       const unsigned long long int config_dump_steps,
                                       const unsigned long long int maximum_steps,
                                       const double initial_temperature,
                                       const std::string &json_coefficients_filename)
    : McAbstract(cfg::GenerateSoluteConfig(factors, solvent_element, solute_atom_count),
                 log_dump_steps,
                 config_dump_steps,
                 maximum_steps,
                 0,
                 0,
                 0,
                 0,
                 initial_temperature,
                 GetElementSetFromSolventAndSolute(solvent_element, solute_atom_count),
                 json_coefficients_filename,
                 "sa_log.txt"),
      energy_change_predictor_(
          json_coefficients_filename, config_, GetElementSetFromSolventAndSolute(solvent_element, solute_atom_count)),
      atom_index_selector_(0, config_.GetNumAtoms() - 1),
      initial_temperature_(initial_temperature),
      reheat_trigger_steps_(std::max<unsigned long long int>(1ULL, static_cast<unsigned long long>(maximum_steps_ * reheat_trigger_ratio_))),
      reheat_cooldown_steps_(std::max<unsigned long long int>(1ULL, static_cast<unsigned long long>(maximum_steps_ * reheat_cooldown_ratio_))){
  ofs_.precision(16);
  const auto energy_predictor = pred::EnergyPredictor(
      json_coefficients_filename, GetElementSetFromSolventAndSolute(solvent_element, solute_atom_count));
  auto chemical_potential = energy_predictor.GetChemicalPotential(solvent_element);

  double energy_change_solution_to_pure_solvent = 0;
  for (const auto &[element, count]: solute_atom_count) {
    energy_change_solution_to_pure_solvent += chemical_potential[element] * static_cast<double>(count);
  }
  const double energy_change_cluster_to_pure_solvent =
      energy_predictor.GetEnergy(config_) - energy_predictor.GetEnergy(cfg::GenerateFCC(factors, solvent_element));
  energy_ = energy_change_cluster_to_pure_solvent - energy_change_solution_to_pure_solvent;
  std::cout << "initial_energy = " << energy_ << std::endl;
  // initialize schedule state related to energy tracking
  lowest_energy_ = energy_;
  recent_best_energy_ = energy_;
  last_improvement_step_ = steps_;
  // stage concept removed; schedule is fully smooth
}

void SimulatedAnnealing::Dump() const {
  if (steps_ == 0) {
    ofs_ << "steps\ttemperature\tenergy\tlowest_energy\tabsolute_energy" << std::endl;
  }

  if (energy_ < lowest_energy_ - kEpsilon) {
    lowest_energy_ = energy_;
    config_.WriteConfig("lowest_energy.cfg.gz");
    ofs_ << steps_ << '\t' << temperature_ << '\t' << energy_ << '\t' << lowest_energy_ << '\t' << absolute_energy_
         << std::endl;
  }

  if (steps_ % log_dump_steps_ == 0) {
    ofs_ << steps_ << '\t' << temperature_ << '\t' << energy_ << '\t' << lowest_energy_ << '\t' << absolute_energy_
         << std::endl;
  }
//  if (steps_ % maximum_steps_ == 0 && steps_ != 0) {
//    config_.WriteConfig(std::to_string(static_cast<int>(temperature_)) + "K.cfg.gz");
//  }
}
void SimulatedAnnealing::UpdateTemperature() {
  // Dynamic target band shrinks with progress: early (acc_early_.first..acc_early_.second)
  // → late (acc_late_.first..acc_late_.second)
  const double progress = static_cast<double>(steps_) / static_cast<double>(std::max<unsigned long long int>(1ULL, maximum_steps_));
  const double acc_max = acc_early_.second + (acc_late_.second - acc_early_.second) * progress;
  const double acc_min = acc_early_.first + (acc_late_.first - acc_early_.first) * progress;

  // Acceptance window based temperature adjustment (thermostat-like)
  if (window_trials_ >= window_size_) {
    const double acc = static_cast<double>(window_accepts_) / static_cast<double>(window_trials_);
    (void)acc_min; // not heating by default; reserved for future use
    if (acc > acc_max) {
      temperature_ *= cool_down_factor_;
    }
    window_trials_ = 0;
    window_accepts_ = 0;
  }

  // Reheat on stagnation: require (no recent improvement) AND (acceptance too low)
  const double acc_est = (window_trials_ > 0u) ? (static_cast<double>(window_accepts_) / static_cast<double>(window_trials_)) : 1.0;
  if (reheats_done_ < max_reheats_ &&
      (steps_ - last_improvement_step_ >= reheat_trigger_steps_) &&
      (steps_ - last_reheat_step_ >= reheat_cooldown_steps_) &&
      (acc_est < acc_min)) {
    temperature_ *= reheat_factor_;
    last_improvement_step_ = steps_;
    last_reheat_step_ = steps_;
    recent_best_energy_ = energy_;
    ++reheats_done_;
  }

  // Continuous baseline cooling to ensure monotone descent trend
  temperature_ *= std::exp(-baseline_c_ / static_cast<double>(std::max<unsigned long long int>(1ULL, maximum_steps_)));

  // update beta
  beta_ = 1.0 / constants::kBoltzmann / std::max(temperature_, 1e-12);
}
std::pair<size_t, size_t> SimulatedAnnealing::GenerateLatticeIdJumpPair() {
  size_t lattice_id1, lattice_id2;
  do {
    lattice_id1 = atom_index_selector_(generator_);
    lattice_id2 = atom_index_selector_(generator_);
  } while (config_.GetElementAtLatticeId(lattice_id1)
      == config_.GetElementAtLatticeId(lattice_id2));
  return {lattice_id1, lattice_id2};
}
bool SimulatedAnnealing::SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                     const double dE) {
  if (dE < 0) {
    config_.LatticeJump(lattice_id_jump_pair);
    energy_ += dE;
    absolute_energy_ += dE;
    return true;
  } else {
    double possibility = std::exp(-dE * beta_);
    double random_number = unit_distribution_(generator_);
    if (random_number < possibility) {
      config_.LatticeJump(lattice_id_jump_pair);
      energy_ += dE;
      absolute_energy_ += dE;
      return true;
    }
  }
  return false;
}
void SimulatedAnnealing::Simulate() {
  while (steps_ <= maximum_steps_) {
    auto lattice_id_jump_pair = GenerateLatticeIdJumpPair();
    auto dE = energy_change_predictor_.GetDeFromLatticeIdPair(config_, lattice_id_jump_pair);

    Dump();

    // perform Metropolis step at current temperature
    const bool accepted = SelectEvent(lattice_id_jump_pair, dE);

    // update acceptance window stats
    ++window_trials_;
    if (accepted) {
      ++window_accepts_;
      // improvement detection relative to recent baseline (not global lowest)
      if (energy_ < recent_best_energy_ - kEpsilon) {
        recent_best_energy_ = energy_;
        last_improvement_step_ = steps_;
      }
    }

    // statistics and logging
    thermodynamic_averaging_.AddEnergy(energy_);

    // update schedule for next step
    UpdateTemperature();

    ++steps_;
  }
}

}    // namespace mc

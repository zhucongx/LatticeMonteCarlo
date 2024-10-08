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
      lowest_energy_(0),
      initial_temperature_(initial_temperature){
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
  static const double alpha = 1. - 3. / static_cast<double>(maximum_steps_);
  temperature_ = initial_temperature_ * std::pow(alpha, steps_);
  beta_ = 1.0 / constants::kBoltzmann / temperature_;
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
void SimulatedAnnealing::SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                     const double dE) {
  if (dE < 0) {
    config_.LatticeJump(lattice_id_jump_pair);
    energy_ += dE;
    absolute_energy_ += dE;
  } else {
    double possibility = std::exp(-dE * beta_);
    double random_number = unit_distribution_(generator_);
    if (random_number < possibility) {
      config_.LatticeJump(lattice_id_jump_pair);
      energy_ += dE;
      absolute_energy_ += dE;
    }
  }
}
void SimulatedAnnealing::Simulate() {
  while (steps_ <= maximum_steps_) {
    auto lattice_id_jump_pair = GenerateLatticeIdJumpPair();
    auto dE = energy_change_predictor_.GetDeFromLatticeIdPair(config_, lattice_id_jump_pair);
    thermodynamic_averaging_.AddEnergy(energy_);
    Dump();
    UpdateTemperature();
    SelectEvent(lattice_id_jump_pair, dE);
    ++steps_;
  }
}

}    // namespace mc
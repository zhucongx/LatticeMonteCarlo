#include "CanonicalMcStepT.h"
#include <utility>
#include <chrono>
#include <omp.h>
#include "EnergyPredictor.h"
namespace mc {
// static std::set<Element> GetElementSetFromSolventAndSolute(
//     Element solvent_element, const std::set<Element> &solute_element_set) {
//
//   std::set<Element> element_set;
//   element_set.insert(solvent_element);
//   for (const auto &solute_element: solute_element_set) {
//     element_set.insert(solute_element);
//   }
//   return element_set;
// }
CanonicalMcStepT::CanonicalMcStepT(cfg::Config config,
                                   const std::set<Element> &element_set,
                                   const unsigned long long int log_dump_steps,
                                   const unsigned long long int config_dump_steps,
                                   const unsigned long long int maximum_steps,
                                   const unsigned long long int thermodynamic_averaging_steps,
                                   double initial_temperature,
                                   double decrement_temperature,
                                   const std::string &json_coefficients_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_steps_(maximum_steps),
      thermodynamic_averaging_(thermodynamic_averaging_steps),
      initial_temperature_(initial_temperature),
      decrement_temperature_(decrement_temperature),
      temperature_(initial_temperature),
      beta_(1 / constants::kBoltzmannConstant / initial_temperature_),
      energy_predictor_(json_coefficients_filename,
                        config_,
                        element_set),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      atom_index_selector_(0, config_.GetNumAtoms() - 1) {
  pred::EnergyPredictor total_energy_predictor(json_coefficients_filename, element_set);
#pragma omp parallel master default(none) shared(std::cout)
  {
    std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
  }
  std::ofstream ofs("cmc_log.txt", std::ofstream::out);
  ofs.precision(16);
  ofs << "initial_energy = " << total_energy_predictor.GetEnergy(config_) << std::endl;
  ofs << "steps\ttemperature\tenergy\taverage_energy\n";
}
std::pair<size_t, size_t> CanonicalMcStepT::GenerateAtomIdJumpPair() {
  size_t atom_id1, atom_id2;
  do {
    atom_id1 = atom_index_selector_(generator_);
    atom_id2 = atom_index_selector_(generator_);
  } while (atom_id1 == atom_id2);
  return {atom_id1, atom_id2};
}

void CanonicalMcStepT::UpdateTemperature() {
  if (steps_ % maximum_steps_ == 0 && steps_ != 0) {
    config_.WriteConfig("end_" + std::to_string(static_cast<int>(temperature_)) + "K.cfg",
                        false);
    temperature_ -= decrement_temperature_;
    beta_ = 1.0 / constants::kBoltzmannConstant / temperature_;
  }
}

void CanonicalMcStepT::Dump(std::ofstream &ofs) {
  unsigned long long int log_dump_steps;
  if (steps_ > 10 * log_dump_steps_) {
    log_dump_steps = log_dump_steps_;
  } else {
    log_dump_steps = static_cast<unsigned long long int>(
        std::pow(10, static_cast<unsigned long long int>(std::log10(
            steps_ + 1) - 1)));
    log_dump_steps = std::max(log_dump_steps, static_cast<unsigned long long int>(1));
    log_dump_steps = std::min(log_dump_steps, log_dump_steps_);
  }
  if (steps_ % log_dump_steps == 0) {
    ofs << steps_ << '\t' << temperature_ << '\t' << energy_
        << thermodynamic_averaging_.GetThermodynamicAverage(temperature_) << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteConfig(std::to_string(steps_) + ".cfg", false);
  }
}

void CanonicalMcStepT::Simulate() {
  std::ofstream ofs("cmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs.precision(16);
  auto t1 = std::chrono::high_resolution_clock::now();
  while (steps_ <= maximum_steps_
      * static_cast<unsigned long long int>(initial_temperature_ / decrement_temperature_ + 1)) {
    thermodynamic_averaging_.AddEnergy(energy_);
    Dump(ofs);
    UpdateTemperature();
    auto atom_id_jump_pair = GenerateAtomIdJumpPair();
    auto dE = energy_predictor_.GetDeFromAtomIdPair(config_, atom_id_jump_pair);
    if (dE < 0) {
      config_.AtomJump(atom_id_jump_pair);
      energy_ += dE;
    } else {
      double possibility = std::exp(-dE * beta_);
      double random_number = one_distribution_(generator_);
      if (random_number < possibility) {
        config_.AtomJump(atom_id_jump_pair);
        energy_ += dE;
      }
    }
    ++steps_;
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Canonical Monte Carlo finished in " << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds.\n";
}

} // mc

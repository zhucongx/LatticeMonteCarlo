#include "CanonicalMcNearStepT.h"
#include <utility>
#include <chrono>
#include <omp.h>
#include "TotalEnergyPredictor.h"
namespace ansys {
constexpr double kBoltzmannConstant = 8.617333262145e-5;
static std::set<Element> GetElementSetFromSolventAndSolute(
    Element solvent_element, const std::set<Element> &solute_element_set) {

  std::set<Element> element_set;
  element_set.insert(solvent_element);
  for (const auto &solute_element: solute_element_set) {
    element_set.insert(solute_element);
  }
  return element_set;
}
CanonicalMcNearStepT::CanonicalMcNearStepT(cfg::Config config,
                                           Element solvent_element,
                                           const std::set<Element> &solute_element_set,
                                           const unsigned long long int log_dump_steps,
                                           const unsigned long long int config_dump_steps,
                                           const unsigned long long int maximum_number,
                                           double initial_temperature,
                                           double decrement_temperature,
                                           const std::string &json_coefficients_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      initial_temperature_(initial_temperature),
      decrement_temperature_(decrement_temperature),
      temperature_(initial_temperature),
      beta_(1 / kBoltzmannConstant / initial_temperature_),
      energy_predictor_(json_coefficients_filename,
                        config_,
                        GetElementSetFromSolventAndSolute(solvent_element, solute_element_set)),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      atom_index_selector_(0, config_.GetNumAtoms() - 1),
      neighbor_index_selector_(0, constants::kNumFirstNearestNeighbors - 1) {
  pred::TotalEnergyPredictor total_energy_predictor(
      json_coefficients_filename,
      GetElementSetFromSolventAndSolute(solvent_element, solute_element_set));
#pragma omp parallel master default(none) shared(std::cout)
  {
    std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
  }
  std::ofstream ofs("cmc_log.txt", std::ofstream::out);
  ofs.precision(16);
  ofs << "initial_energy = " << total_energy_predictor.GetEnergy(config_) << std::endl;
  ofs << "steps\tenergy\ttemperature\n";

}
std::pair<size_t, size_t> CanonicalMcNearStepT::GenerateAtomIdJumpPair() {
  size_t atom_id1 = atom_index_selector_(generator_);
  size_t lattice_id1 = config_.GetLatticeIdFromAtomId(atom_id1);
  size_t lattice_id2 = config_.GetFirstNeighborsAdjacencyList().at(
      lattice_id1).at(neighbor_index_selector_(generator_));
  size_t atom_id2 = config_.GetAtomIdFromLatticeId(lattice_id2);
  return {atom_id1, atom_id2};
}

void CanonicalMcNearStepT::UpdateTemperature() {
  if (steps_ % maximum_number_ == 0 && steps_ != 0) {
    config_.WriteCfg("end_" + std::to_string(static_cast<int>(temperature_)) + "K.cfg",
                     false);
    temperature_ -= decrement_temperature_;
    beta_ = 1.0 / kBoltzmannConstant / temperature_;
  }
}

void CanonicalMcNearStepT::Dump(std::ofstream &ofs) {
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
    ofs << steps_ << '\t' << energy_ << '\t' << temperature_ << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteCfg(std::to_string(steps_) + ".cfg", false);
  }
}

void CanonicalMcNearStepT::Simulate() {
  std::ofstream ofs("cmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs.precision(16);
  auto t1 = std::chrono::high_resolution_clock::now();
  while (steps_ <= maximum_number_
      * static_cast<unsigned long long int>(initial_temperature_ / decrement_temperature_ + 1)) {
    Dump(ofs);
    UpdateTemperature();
    auto atom_id_jump_pair = GenerateAtomIdJumpPair();
    auto dE = energy_predictor_.GetDiffFromAtomIdPair(config_, atom_id_jump_pair);
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

  config_.WriteCfg("final.cfg", false);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Canonical Monte Carlo finished in " << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds.\n";
}

} // namespace ansys

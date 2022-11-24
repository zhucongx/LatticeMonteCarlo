#include "CanonicalMc.h"
#include <utility>
#include <chrono>
#include "EnergyPredictor.h"
namespace mc {

CanonicalMc::CanonicalMc(cfg::Config config,
                         const std::set<Element> &element_set,
                         const unsigned long long int log_dump_steps,
                         const unsigned long long int config_dump_steps,
                         const unsigned long long int maximum_number,
                         double temperature,
                         const std::string &json_coefficients_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      beta_(1 / temperature / kBoltzmannConstant),
      energy_predictor_(json_coefficients_filename,
                        config_,
                        element_set),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      atom_index_selector_(0, config_.GetNumAtoms() - 1),
      neighbor_index_selector_(0, constants::kNumFirstNearestNeighbors - 1) {
  pred::EnergyPredictor total_energy_predictor(json_coefficients_filename, element_set);

  std::ofstream ofs("cmc_log.txt", std::ofstream::out);
  ofs.precision(16);
  ofs << "initial_energy = " << total_energy_predictor.GetEnergy(config_) << " T = "
      << temperature << std::endl;
  ofs << "steps\tenergy\n";

}
std::pair<size_t, size_t> CanonicalMc::GenerateAtomIdJumpPair() {
  size_t atom_id1 = atom_index_selector_(generator_);
  size_t lattice_id1 = config_.GetLatticeIdFromAtomId(atom_id1);
  size_t lattice_id2 = config_.GetFirstNeighborsAdjacencyList().at(
      lattice_id1).at(neighbor_index_selector_(generator_));
  size_t atom_id2 = config_.GetAtomIdFromLatticeId(lattice_id2);
  return {atom_id1, atom_id2};
}
void CanonicalMc::Dump(std::ofstream &ofs) {
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
    ofs << steps_ << '\t' << energy_ << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteConfig(std::to_string(steps_) + ".cfg", false);
  }
}

void CanonicalMc::Simulate() {
  std::ofstream ofs("cmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs.precision(8);
  auto t1 = std::chrono::high_resolution_clock::now();
  while (steps_ <= maximum_number_) {
    Dump(ofs);
    auto atom_id_jump_pair = GenerateAtomIdJumpPair();
    auto dE = energy_predictor_.GetDeFromAtomIdPair(
        config_, atom_id_jump_pair);
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

  config_.WriteConfig("final.cfg", false);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Canonical Monte Carlo finished in " << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds.\n";
}

} // mc

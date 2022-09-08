#include "CanonicalMC.h"
#include <utility>
#include <chrono>
#include "TotalEnergyPredictor.h"
namespace ansys {
constexpr double kBoltzmannConstant = 8.617333262145e-5;

CanonicalMC::CanonicalMC(cfg::Config config,
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
  pred::TotalEnergyPredictor total_energy_predictor(json_coefficients_filename,
                                                    element_set);

  std::ofstream ofs("cmc_log.txt", std::ofstream::out);
  ofs << "initial_energy = " << total_energy_predictor.GetEnergy(config_) << ", T = "
      << temperature << std::endl;
  ofs << "steps\tenergy\n";

}
std::pair<size_t, size_t> CanonicalMC::GenerateAtomIdJumpPair() {
  size_t atom_id1 = atom_index_selector_(generator_);
  size_t lattice_id1 = config_.GetLatticeIdFromAtomId(atom_id1);
  size_t lattice_id2 = config_.GetFirstNeighborsAdjacencyList().at(
      lattice_id1).at(neighbor_index_selector_(generator_));
  size_t atom_id2 = config_.GetAtomIdFromLatticeId(lattice_id2);
  return {atom_id1, atom_id2};
}
void CanonicalMC::Dump(std::ofstream &ofs) {
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << energy_ << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteCfg(std::to_string(steps_) + ".cfg", false);
  }
}

void CanonicalMC::Simulate() {
  std::ofstream ofs("cmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs.precision(8);
  auto t1 = std::chrono::high_resolution_clock::now();
  while (steps_ < maximum_number_) {
    Dump(ofs);
    auto atom_id_jump_pair = GenerateAtomIdJumpPair();
    auto dE = energy_predictor_.GetDiffFromAtomIdPair(
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

  config_.WriteCfg("final.cfg", false);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Canonical Monte Carlo finished in " << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds.\n";
}

} // namespace ansys

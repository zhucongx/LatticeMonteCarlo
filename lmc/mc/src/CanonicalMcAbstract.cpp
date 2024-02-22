#include "CanonicalMcAbstract.h"
#include <utility>
#include <chrono>
#include <omp.h>
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
CanonicalMcAbstract::CanonicalMcAbstract(cfg::Config config,
                                         const unsigned long long int log_dump_steps,
                                         const unsigned long long int config_dump_steps,
                                         const unsigned long long int maximum_steps,
                                         const unsigned long long int thermodynamic_averaging_steps,
                                         const unsigned long long int restart_steps,
                                         const double restart_energy,
                                         const double temperature,
                                         const std::set<Element> &element_set,
                                         const std::string &json_coefficients_filename)
    : McAbstract(std::move(config),
                 log_dump_steps,
                 config_dump_steps,
                 maximum_steps,
                 thermodynamic_averaging_steps,
                 restart_steps,
                 restart_energy,
                 0,
                 temperature,
                 element_set,
                 json_coefficients_filename,
                 "cmc_log.txt"),
      energy_change_predictor_(json_coefficients_filename,
                               config_,
                               element_set),
      atom_index_selector_(0, config_.GetNumAtoms() - 1) {
}
std::pair<size_t, size_t> CanonicalMcAbstract::GenerateLatticeIdJumpPair() {
  size_t lattice_id1, lattice_id2;
  do {
    lattice_id1 = atom_index_selector_(generator_);
    lattice_id2 = atom_index_selector_(generator_);
  } while (config_.GetElementAtLatticeId(lattice_id1)
      == config_.GetElementAtLatticeId(lattice_id2));
  return {lattice_id1, lattice_id2};
}

void CanonicalMcAbstract::Dump() const {
  if (is_restarted_) {
    is_restarted_ = false;
    return;
  }
  if (steps_ == 0) {
    config_.WriteLattice("lattice.txt");
    config_.WriteElement("element.txt");
    ofs_ << "steps\ttemperature\tenergy\taverage_energy\tabsolute_energy" << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    // config_.WriteMap("map" + std::to_string(steps_) + ".txt");
    config_.WriteConfig(std::to_string(steps_) + ".cfg.gz");
  }
  if (steps_ == maximum_steps_) {
    config_.WriteConfig("end.cfg");
  }
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
    ofs_ << steps_ << '\t' << temperature_ << '\t' << energy_ << '\t'
         << thermodynamic_averaging_.GetThermodynamicAverage(beta_) << '\t'
         << absolute_energy_ << std::endl;
  }
}
void CanonicalMcAbstract::SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair,
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
} // mc
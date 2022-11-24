#include "SemiGrandCanonicalMcStepT.h"
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
SemiGrandCanonicalMcStepT::SemiGrandCanonicalMcStepT(cfg::Config config,
                                                     const std::set<Element> &element_set,
                                                     const unsigned long long int log_dump_steps,
                                                     const unsigned long long int config_dump_steps,
                                                     const unsigned long long int maximum_steps,
                                                     const unsigned long long int thermodynamic_averaging_steps,
                                                     double initial_temperature,
                                                     double decrement_temperature,
                                                     const std::string &json_coefficients_filename)
    : config_(std::move(config)),
      element_vector_(element_set.begin(), element_set.end()),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_steps_(maximum_steps),
      initial_temperature_(initial_temperature),
      decrement_temperature_(decrement_temperature),
      temperature_(initial_temperature),
      beta_(1 / constants::kBoltzmannConstant / initial_temperature_),
      energy_predictor_(json_coefficients_filename,
                        config_,
                        element_set),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      atom_index_selector_(0, config_.GetNumAtoms() - 1),
      element_index_selector_(0, element_vector_.size() - 1) {
  pred::EnergyPredictor total_energy_predictor(json_coefficients_filename, element_set);
  chemical_potential_ = total_energy_predictor.GetChemicalPotential();
#pragma omp parallel master default(none) shared(std::cout)
  {
    std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
  }
  std::ofstream ofs("sgcmc_log.txt", std::ofstream::out);
  ofs.precision(16);
  double total_energy = total_energy_predictor.GetEnergy(config_);
  for (const auto &element_vector: config_.GetElementAtomIdVectorMap()) {
    total_energy -= chemical_potential_[element_vector.first]
        * static_cast<double>(element_vector.second.size());
  }

  ofs << "initial_energy = " << total_energy << std::endl;
  ofs << "steps\ttemperature\tenergy\taverage_energy\n";
}
std::pair<size_t, Element> SemiGrandCanonicalMcStepT::GenerateAtomIdChangeSite() {
  size_t atom_id;
  Element element;
  do {
    atom_id = atom_index_selector_(generator_);
    element = element_vector_[element_index_selector_(generator_)];
  } while (config_.GetElementAtAtomId(atom_id) == ElementName::X);
  return {atom_id, element};
}

void SemiGrandCanonicalMcStepT::UpdateTemperature() {
  if (steps_ % maximum_steps_ == 0 && steps_ != 0) {
    config_.WriteConfig("end_" + std::to_string(static_cast<int>(temperature_)) + "K.cfg",
                        false);
    temperature_ -= decrement_temperature_;
    beta_ = 1.0 / constants::kBoltzmannConstant / temperature_;
  }
}

void SemiGrandCanonicalMcStepT::Dump(std::ofstream &ofs) {
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
    ofs << steps_ << '\t' << temperature_ << '\t' << energy_ << '\t'  << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteConfig(std::to_string(steps_) + ".cfg", false);
  }
}

void SemiGrandCanonicalMcStepT::Simulate() {
  std::ofstream ofs("sgcmc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs.precision(16);
  auto t1 = std::chrono::high_resolution_clock::now();
  while (steps_ <= maximum_steps_
      * static_cast<unsigned long long int>(initial_temperature_ / decrement_temperature_ + 1)) {
    Dump(ofs);
    UpdateTemperature();
    auto [atom_id, new_element] = GenerateAtomIdChangeSite();
    auto old_element = config_.GetElementAtAtomId(atom_id);
    const auto dE = energy_predictor_.GetDeFromAtomIdSite(config_, atom_id, new_element)
        - (chemical_potential_[new_element] - chemical_potential_[old_element]);
    // std::cerr << "element = " << new_element.GetString() << " dE = " << dE << std::endl;
    if (dE < 0) {
      config_.ChangeAtomElementTypeAtAtom(atom_id, new_element);
      energy_ += dE;
    } else {
      double possibility = std::exp(-dE * beta_);
      double random_number = one_distribution_(generator_);
      if (random_number < possibility) {
        config_.ChangeAtomElementTypeAtAtom(atom_id, new_element);
        energy_ += dE;
      }
    }
    ++steps_;
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Semi Grand Canonical Monte Carlo finished in " << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds.\n";
}

} // mc
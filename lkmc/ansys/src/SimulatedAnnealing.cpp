#include "SimulatedAnnealing.h"
#include <utility>
#include <chrono>
#include "TotalEnergyPredictor.h"
namespace ansys {
static std::set<Element> GetElementSetFromSolventAndSolute(
    Element solvent_element,
    const std::map<Element, size_t> &solute_atom_count) {
  std::set<Element> element_hashset;
  element_hashset.insert(solvent_element);
  for (const auto &solute_element_count: solute_atom_count) {
    element_hashset.insert(solute_element_count.first);
  }
  return element_hashset;
}
static std::vector<size_t> GetSoluteAtomIdVector(
    const cfg::Config &config,
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
                                       const unsigned long long int maximum_number,
                                       unsigned long long int early_stop_number,
                                       const std::string &json_coefficients_filename)
    : config_(cfg::GenerateSoluteConfig(factors, solvent_element, solute_atom_count)),
      solute_atom_id_vector_(GetSoluteAtomIdVector(config_, solute_atom_count)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      early_stop_number_(early_stop_number),
      steps_(0),
      count_(0),
      energy_predictor_(json_coefficients_filename,
                        config_,
                        GetElementSetFromSolventAndSolute(solvent_element, solute_atom_count)),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      solute_atom_selector_(0, solute_atom_id_vector_.size() - 1),
      neighbor_index_selector_(0, constants::kNumFirstNearestNeighbors - 1) {
  // pred::TotalEnergyPredictor total_energy_predictor(json_coefficients_filename,
  //                                                   GetElementSetFromSolventAndSolute(
  //                                                       solvent_element,
  //                                                       solute_atom_count));
  // energy_ = total_energy_predictor.GetEnergy(config_);
  // lowest_energy_ = energy_;
}
std::pair<size_t, size_t> SimulatedAnnealing::GenerateAtomIdJumpPair() {
  size_t atom_id1 =solute_atom_id_vector_.at(solute_atom_selector_(generator_));
  size_t lattice_id1 = config_.GetLatticeIdFromAtomId(atom_id1);
  size_t lattice_id2 = config_.GetFirstNeighborsAdjacencyList().at(
      lattice_id1).at(neighbor_index_selector_(generator_));
  size_t atom_id2 = config_.GetAtomIdFromLatticeId(lattice_id2);
  return {atom_id1, atom_id2};
}
void SimulatedAnnealing::Dump(std::ofstream &ofs) {
  if (energy_ < lowest_energy_) {
    lowest_energy_ = energy_;
    config_.WriteCfg("lowest_energy.cfg", false);
    ofs << steps_ << '\t' << energy_ << '\t' << lowest_energy_ << '\t'
        << temperature_ / kBoltzmannConstant << std::endl;
  }
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << energy_ << '\t' << lowest_energy_ << '\t'
        << temperature_ / kBoltzmannConstant << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteCfg(std::to_string(steps_) + ".cfg", false);
  }
}

void SimulatedAnnealing::Simulate() {
  std::ofstream ofs("mc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs << "steps\tenergy\tlowest_energy\ttemperature\n";
  ofs.precision(8);

  while (steps_ < maximum_number_) {
    temperature_ = initial_temperature_ / std::log(2 + steps_);

    auto atom_id_jump_pair = GenerateAtomIdJumpPair();
    auto dE = energy_predictor_.GetDiffFromAtomIdPair(
        config_, atom_id_jump_pair);
    if (dE < 0) {
      config_.AtomJump(atom_id_jump_pair);
      energy_ += dE;
    } else {
      double possibility = std::exp(-dE / temperature_);
      double random_number = one_distribution_(generator_);
      if (random_number < possibility) {
        config_.AtomJump(atom_id_jump_pair);
        energy_ += dE;
      }
    }
    // if (EarlyStop()) {
    //   Dump(ofs);
    //   break;
    // }
    Dump(ofs);
    ++steps_;
  }
}
double SimulatedAnnealing::MoveOneStep() {
  auto atom_id_jump_pair = GenerateAtomIdJumpPair();
  auto dE = energy_predictor_.GetDiffFromAtomIdPair(
      config_, atom_id_jump_pair);
  config_.AtomJump(atom_id_jump_pair);
  return dE;
}
bool SimulatedAnnealing::EarlyStop() {
  if (energy_ >= lowest_energy_) {
    ++count_;
  } else {
    count_ = 0;
  }
  if (count_ > early_stop_number_) {
    return true;
  }
  return false;
}

} // namespace ansys

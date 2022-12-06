#include "CanonicalMcOmp.h"
#include <omp.h>
namespace mc {
CanonicalMcOmp::CanonicalMcOmp(cfg::Config config,
                               const unsigned long long int log_dump_steps,
                               const unsigned long long int config_dump_steps,
                               const unsigned long long int maximum_steps,
                               const unsigned long long int thermodynamic_averaging_steps,
                               const double initial_temperature,
                               const double decrement_temperature,
                               const std::set<Element> &element_set,
                               const std::string &json_coefficients_filename)
    : CanonicalMcAbstract(std::move(config),
                          log_dump_steps,
                          config_dump_steps,
                          maximum_steps,
                          thermodynamic_averaging_steps,
                          initial_temperature,
                          decrement_temperature,
                          element_set,
                          json_coefficients_filename) {
  if (world_size_ != 1) {
    std::cout << "Must use 1 precesses. Terminating...\n" << std::endl;
    MPI_Finalize();
    exit(0);
  }
#pragma omp parallel master default(none) shared(std::cout)
  {
    std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
    event_vector_.resize(static_cast<size_t>(omp_get_num_threads()));
  }
}
void CanonicalMcOmp::BuildEventVector() {
  std::unordered_set<size_t> unavailable_position{};
  std::pair<size_t, size_t> lattice_id_jump_pair;
  for (auto &event: event_vector_) {
    do {
      lattice_id_jump_pair = GenerateLatticeIdJumpPair();
    } while (unavailable_position.find(lattice_id_jump_pair.first) != unavailable_position.end() ||
        unavailable_position.find(lattice_id_jump_pair.second) != unavailable_position.end());
    for (auto selected_lattice_index: {lattice_id_jump_pair.first, lattice_id_jump_pair.second}) {
      unavailable_position.emplace(selected_lattice_index);
      std::copy(config_.GetFirstNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config_.GetFirstNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position,
                              unavailable_position.end()));
      std::copy(config_.GetSecondNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config_.GetSecondNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position,
                              unavailable_position.end()));
      std::copy(config_.GetThirdNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config_.GetThirdNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position,
                              unavailable_position.end()));
    }
    event = {lattice_id_jump_pair, 0};
  }
#pragma omp parallel for default(none)
  for (auto &event: event_vector_) {
    event.second = energy_change_predictor_.GetDeFromLatticeIdPair(config_, event.first);
  }
}
void CanonicalMcOmp::Simulate() {
  while (steps_ <= maximum_steps_
      * static_cast<unsigned long long int>(initial_temperature_ / decrement_temperature_ + 1)) {
    BuildEventVector();
    for (auto [lattice_id_jump_pair, dE]: event_vector_) {
      thermodynamic_averaging_.AddEnergy(energy_);
      Dump();
      UpdateTemperature();
      SelectEvent(lattice_id_jump_pair, dE);
      ++steps_;
    }
  }
}
} // mc
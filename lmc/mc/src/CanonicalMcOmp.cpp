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
    num_threads_ = static_cast<size_t>(omp_get_num_threads());
  }
  event_vector_.reserve(static_cast<size_t>(omp_get_num_threads()));
}
void CanonicalMcOmp::BuildEventVector() {
  event_vector_.clear();
  unavailable_position_.clear();

  std::pair<size_t, size_t> lattice_id_jump_pair;
  for (size_t i = 0; i < num_threads_; ++i) {
    size_t ct = 0;
    do {
      lattice_id_jump_pair = GenerateLatticeIdJumpPair();
      ct++;
      if (ct == 50) {
        break;
      }
    } while (unavailable_position_.find(lattice_id_jump_pair.first) != unavailable_position_.end()
        || unavailable_position_.find(lattice_id_jump_pair.second) != unavailable_position_.end());
    if (ct == 50) {
      break;
    }
    for (auto selected_lattice_index: {lattice_id_jump_pair.first, lattice_id_jump_pair.second}) {
      unavailable_position_.emplace(selected_lattice_index);
      std::copy(config_.GetFirstNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config_.GetFirstNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position_,
                              unavailable_position_.end()));
      std::copy(config_.GetSecondNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config_.GetSecondNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position_,
                              unavailable_position_.end()));
      std::copy(config_.GetThirdNeighborsAdjacencyList().at(selected_lattice_index).begin(),
                config_.GetThirdNeighborsAdjacencyList().at(selected_lattice_index).end(),
                std::inserter(unavailable_position_,
                              unavailable_position_.end()));
    }
    event_vector_.emplace_back(lattice_id_jump_pair, 0);
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
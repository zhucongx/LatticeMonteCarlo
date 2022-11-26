#include "KineticMcFirstOmp.h"
#include <utility>
#include <chrono>
namespace mc {

KineticMcFirstOmp::KineticMcFirstOmp(cfg::Config config,
                                     unsigned long long int log_dump_steps,
                                     unsigned long long int config_dump_steps,
                                     unsigned long long int maximum_steps,
                                     unsigned long long int thermodynamic_averaging_steps,
                                     double temperature,
                                     const std::set<Element> &element_set,
                                     unsigned long long int restart_steps,
                                     double restart_energy,
                                     double restart_time,
                                     const std::string &json_coefficients_filename)
    : KineticMcFirstAbstract(std::move(config),
                             log_dump_steps,
                             config_dump_steps,
                             maximum_steps,
                             thermodynamic_averaging_steps,
                             temperature,
                             element_set,
                             restart_steps,
                             restart_energy,
                             restart_time,
                             json_coefficients_filename) {
  if (world_size_ != 1) {
    std::cout << "Must use " << kEventListSize << " precesses. Terminating...\n" << std::endl;
    MPI_Finalize();
    exit(0);
  }
#pragma omp parallel master default(none) shared(std::cout)
  {
    std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
  }
}
KineticMcFirstOmp::~KineticMcFirstOmp() = default;
void KineticMcFirstOmp::BuildEventList() {
  total_rate_k_ = 0.0;
  const auto neighbor_vacancy_id_vector =
      config_.GetFirstNeighborsAdjacencyList()[vacancy_lattice_id_];
#pragma omp parallel for default(none) shared(vacancy_lattice_id_, neighbor_vacancy_id_vector) reduction(+: total_rate_k_)
  for (size_t i = 0; i < kEventListSize; ++i) {
    const auto neighbor_vacancy_id = neighbor_vacancy_id_vector[i];
    JumpEvent lattice_jump_event(
        {vacancy_lattice_id_, neighbor_vacancy_id},
        energy_predictor_.GetBarrierAndDiffFromLatticeIdPair(
            config_, {vacancy_lattice_id_, neighbor_vacancy_id}),
        beta_);
    total_rate_k_ += lattice_jump_event.GetForwardRate();
    event_k_i_list_[i] = std::move(lattice_jump_event);
  }
  double cumulative_probability = 0.0;
  for (auto &event_it: event_k_i_list_) {
    event_it.CalculateProbability(total_rate_k_);
    cumulative_probability += event_it.GetProbability();
    event_it.SetCumulativeProbability(cumulative_probability);
  }
}
double KineticMcFirstOmp::CalculateTime() {
  static std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
  return -std::log(uniform_distribution(generator_)) / total_rate_k_ / constants::kPrefactor;
}
} // mc
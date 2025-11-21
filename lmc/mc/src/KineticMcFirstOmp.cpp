#include "KineticMcFirstOmp.h"
#include <utility>
#include <chrono>

namespace mc {

KineticMcFirstOmp::KineticMcFirstOmp(cfg::Config config,
                                     const unsigned long long int log_dump_steps,
                                     const unsigned long long int config_dump_steps,
                                     const unsigned long long int maximum_steps,
                                     const unsigned long long int thermodynamic_averaging_steps,
                                     const unsigned long long int restart_steps,
                                     const double restart_energy,
                                     const double restart_time,
                                     const double temperature,
                                     const std::set<Element> &element_set,
                                     const std::string &json_coefficients_filename,
                                     const std::string &time_temperature_filename,
                                     const bool is_rate_corrector,
                                     const bool is_early_stop,
                                     const bool is_solute_disp)
    : KineticMcFirstAbstract(std::move(config),
                             log_dump_steps,
                             config_dump_steps,
                             maximum_steps,
                             thermodynamic_averaging_steps,
                             restart_steps,
                             restart_energy,
                             restart_time,
                             temperature,
                             element_set,
                             json_coefficients_filename,
                             time_temperature_filename,
                             is_rate_corrector,
                             is_early_stop,
                             is_solute_disp) {
  if (world_size_ != 1) {
    if (world_rank_ == 0) {
      std::cout << "Must use 1 process for KineticMcFirstOmp. Terminating...\n" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#pragma omp parallel default(none) shared(std::cout)
  {
#pragma omp master
    {
      std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
    }
  }
}
KineticMcFirstOmp::~KineticMcFirstOmp() = default;
  // TODO(perf): This outer OMP over 12 neighbors often contends with the LRU cache
  // locks inside the predictor. Prefer a sequential loop here and keep the predictor's
  // internal parallel sections (dE/D/Ks) enabled, or use thread-local LRU caches.
  // Consider also reusing the previous step's reverse event to avoid one prediction.
void KineticMcFirstOmp::BuildEventList() {
  total_rate_k_ = 0.0;
  const auto neighbor_lattice_id_vector = config_.GetFirstNeighborsAdjacencyList()[vacancy_lattice_id_];
#pragma omp parallel default(none) shared(neighbor_lattice_id_vector) reduction(+: total_rate_k_)
  {
#pragma omp for
    for (size_t i = 0; i < kEventListSize; ++i) {
      const auto neighbor_lattice_id = neighbor_lattice_id_vector[i];
      JumpEvent lattice_jump_event(
          {vacancy_lattice_id_, neighbor_lattice_id},
          vacancy_migration_predictor_lru_.GetBarrierAndDiffFromLatticeIdPair(
              config_, {vacancy_lattice_id_, neighbor_lattice_id}),
          beta_);
      total_rate_k_ += lattice_jump_event.GetForwardRate();
      event_k_i_list_[i] = std::move(lattice_jump_event);
    }
  }
  double cumulative_probability = 0.0;
  for (auto &event_it : event_k_i_list_) {
    // TODO(perf): Consider sampling over a uniform in [0,total_rate_k_] and
    // linearly scanning events by accumulating forward rates, so we can skip
    // building/using normalized probabilities and cumulative sums here.
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

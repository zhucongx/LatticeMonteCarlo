#include "KineticMcFirstMpi.h"
#include <utility>
#include <chrono>
namespace mc {
KineticMcFirstMpi::KineticMcFirstMpi(cfg::Config config,
                                     unsigned long long int log_dump_steps,
                                     unsigned long long int config_dump_steps,
                                     unsigned long long int maximum_steps,
                                     const unsigned long long int thermodynamic_averaging_steps,
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
  if (world_size_ != kEventListSize) {
    std::cout << "Must use " << kEventListSize << " precesses. Terminating...\n" << std::endl;
    MPI_Finalize();
    exit(0);
  }
  if (world_rank_ == 0) {
    std::cout << "Using " << world_size_ << " processes." << std::endl;
  }
}
KineticMcFirstMpi::~KineticMcFirstMpi() = default;
void KineticMcFirstMpi::BuildEventList() {
  total_rate_k_ = 0;
  const auto neighbor_vacancy_id =
      config_.GetFirstNeighborsAdjacencyList()[vacancy_lattice_id_][static_cast<size_t>(world_rank_)];

  event_k_i_ = JumpEvent(
      {vacancy_lattice_id_, neighbor_vacancy_id},
      energy_predictor_.GetBarrierAndDiffFromAtomIdPair(
          config_, {vacancy_lattice_id_, neighbor_vacancy_id}),
      beta_);
  const double this_rate = event_k_i_.GetForwardRate();
  MPI_Allreduce(&this_rate, &total_rate_k_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  event_k_i_.CalculateProbability(total_rate_k_);

  MPI_Allgather(&event_k_i_,
                sizeof(JumpEvent),
                MPI_BYTE,
                event_k_i_list_.data(),
                sizeof(JumpEvent),
                MPI_BYTE,
                MPI_COMM_WORLD);
  double cumulative_probability = 0.0;
  for (auto &event_it: event_k_i_list_) {
    cumulative_probability += event_it.GetProbability();
    event_it.SetCumulativeProbability(cumulative_probability);
  }
}
double KineticMcFirstMpi::CalculateTime() {
  static std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
  return -std::log(uniform_distribution(generator_)) / total_rate_k_ / constants::kPrefactor;
}
} // mc
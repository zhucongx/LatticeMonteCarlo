#include "KineticMcSecondSerial.h"
namespace mc {
//  j -> k -> i ->l
//       |
// current position
KineticMcSecondSerial::KineticMcSecondSerial(cfg::Config config,
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
                                             const bool is_early_stop,
                                             const bool is_solute_disp,
                                             const unsigned long long int seed)
    : KineticMcSecondAbstract(std::move(config),
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
                              is_early_stop,
                              is_solute_disp,
                              seed) {
  if (world_size_ != 1) {
    if (world_rank_ == 0) {
      std::cout << "Must use 1 process for KineticMcSecondSerial. Terminating...\n" << std::endl;
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

KineticMcSecondSerial::~KineticMcSecondSerial() = default;

void KineticMcSecondSerial::BuildEventList() {
  const auto k_lattice_id = vacancy_lattice_id_;
  const auto &i_lattice_id_list = config_.GetFirstNeighborsAdjacencyList()[k_lattice_id];
  for (size_t ii = 0; ii < kEventListSize; ++ii) {
    const auto i_lattice_id = i_lattice_id_list[ii];
    double total_rate_i = 0.0;
    // The jumps i -> l are those of a vacancy on site i, so the vacancy is moved there and back.
    config_.LatticeJump({k_lattice_id, i_lattice_id});
    for (const auto l_lattice_id: config_.GetFirstNeighborsAdjacencyList()[i_lattice_id]) {
      JumpEvent event_i_l({i_lattice_id, l_lattice_id},
                          vacancy_migration_predictor_lru_.GetBarrierAndDiffFromLatticeIdPair(
                              config_, {i_lattice_id, l_lattice_id}),
                          beta_);
      if (l_lattice_id == k_lattice_id) {
        event_k_i_list_[ii] = event_i_l.GetReverseJumpEvent();
      }
      total_rate_i += event_i_l.GetForwardRate();
    }
    config_.LatticeJump({i_lattice_id, k_lattice_id});
    total_rate_i_list_[ii] = total_rate_i;
  }

  total_rate_k_ = 0.0;
  for (const auto &event_k_i: event_k_i_list_) {
    total_rate_k_ += event_k_i.GetForwardRate();
  }
  for (auto &event_k_i: event_k_i_list_) {
    event_k_i.CalculateProbability(total_rate_k_);
  }
}

double KineticMcSecondSerial::CalculateTime() {
  std::array<SecondOrderSums, kEventListSize> summands_list{};
  SecondOrderSums sums{};
  for (size_t ii = 0; ii < kEventListSize; ++ii) {
    summands_list[ii] = GetSecondOrderSummands(event_k_i_list_[ii], total_rate_i_list_[ii]);
    sums += summands_list[ii];
  }

  // The summands need the first-order probabilities of all twelve events, which the probabilities of the step replace.
  double cumulative_probability = 0.0;
  for (size_t ii = 0; ii < kEventListSize; ++ii) {
    auto &event_k_i = event_k_i_list_[ii];
    event_k_i.SetProbability(GetSecondOrderProbability(event_k_i, summands_list[ii], sums));
    cumulative_probability += event_k_i.GetProbability();
    event_k_i.SetCumulativeProbability(cumulative_probability);
  }
  return GetSecondOrderTime(sums);
}
} // mc

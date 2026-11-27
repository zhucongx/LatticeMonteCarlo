#include "KineticMcSecondOmpi.h"
namespace mc {
//  j -> k -> i ->l
//       |
// current position
KineticMcSecondOmpi::KineticMcSecondOmpi(cfg::Config config,
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
    : KineticMcSecondMpi(std::move(config),
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
  if (world_rank_ == 0) {
#pragma omp parallel default(none) shared(std::cout)
    {
#pragma omp master
      {
        std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
      }
    }
  }
}

void KineticMcSecondOmpi::EvaluateSecondJumps(const size_t k_lattice_id, const size_t i_lattice_id) {
  const auto &l_lattice_id_list = config_.GetFirstNeighborsAdjacencyList()[i_lattice_id];
  total_rate_i_ = 0.0;
#pragma omp parallel default(none) shared(i_lattice_id, k_lattice_id, l_lattice_id_list) reduction(+:total_rate_i_)
  {
#pragma omp for
    for (size_t ii = 0; ii < kEventListSize; ++ii) {
      const auto l_lattice_id = l_lattice_id_list[ii];
      JumpEvent event_i_l({i_lattice_id, l_lattice_id},
                          vacancy_migration_predictor_lru_.GetBarrierAndDiffFromLatticeIdPair(
                              config_, {i_lattice_id, l_lattice_id}),
                          beta_);
      if (l_lattice_id == k_lattice_id) {
        event_k_i_ = event_i_l.GetReverseJumpEvent();
      }
      auto r_i_l = event_i_l.GetForwardRate();
      total_rate_i_ += r_i_l;
    }
  }
}
} // mc

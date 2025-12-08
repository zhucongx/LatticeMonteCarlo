#include "McAbstract.h"
#include "EnergyPredictor.h"
#include <chrono>
#include <mpi.h>
#include <utility>

namespace mc {
McAbstract::McAbstract(cfg::Config config,
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
                       const std::string &log_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_steps_(maximum_steps),
      steps_(restart_steps),
      energy_(restart_energy),
      absolute_energy_(pred::EnergyPredictor(json_coefficients_filename, element_set).GetEnergy(config_)),
      time_(restart_time),
      temperature_(temperature),
      beta_(1.0 / constants::kBoltzmann / temperature_),
      is_restarted_(steps_ > 0),
      thermodynamic_averaging_(thermodynamic_averaging_steps),
      generator_(static_cast<unsigned long long int>(std::chrono::system_clock::now().time_since_epoch().count())),
      unit_distribution_(0.0, 1.0),
      ofs_(log_filename, is_restarted_ ? std::ofstream::app : std::ofstream::out)
{
  // TODO(perf): RNG is on the hot path. Consider replacing mt19937_64 + uniform_real_distribution
  // with a faster generator (pcg/xoshiro) and a lightweight uniform_01 helper to shave per-step cost.
  ofs_.precision(16);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
}
McAbstract::~McAbstract() = default;

}    // namespace mc

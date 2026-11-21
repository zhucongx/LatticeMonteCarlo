#include "McAbstract.h"
#include "EnergyPredictor.h"
#include <algorithm>
#include <cstdint>
#include <mpi.h>
#include <utility>

namespace mc {
namespace {
// A restarted run draws a fresh stream instead of replaying the numbers its previous segment already used.
std::mt19937_64 MakeGenerator(const unsigned long long int seed, const unsigned long long int restart_steps) {
  std::seed_seq sequence{static_cast<std::uint32_t>(seed),
                         static_cast<std::uint32_t>(seed >> 32),
                         static_cast<std::uint32_t>(restart_steps),
                         static_cast<std::uint32_t>(restart_steps >> 32)};
  return std::mt19937_64(sequence);
}

// The energy model only has terms for the vacancy and the species of element_set.
cfg::Config RequireElementsOfModel(cfg::Config config, const std::set<Element> &element_set) {
  for (const auto &element: config.GetElementSetWithoutVacancy()) {
    if (element_set.find(element) == element_set.end()) {
      throw std::runtime_error("The configuration contains " + element.GetString() +
                               ", which is not in element_set");
    }
  }
  return config;
}
}    // namespace

unsigned long long int GetLogDumpInterval(const unsigned long long int steps,
                                          const unsigned long long int log_dump_steps) {
  if (steps > 10 * log_dump_steps) {
    return log_dump_steps;
  }
  // 10^(number of decimal digits of (steps + 1) - 2), i.e. 1 below 99 steps, 10 below 999 steps, and so on.
  unsigned long long int interval = 1;
  for (auto remaining = (steps + 1) / 100; remaining > 0; remaining /= 10) {
    interval *= 10;
  }
  return std::min(interval, log_dump_steps);
}

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
                       const std::string &log_filename,
                       const unsigned long long int seed)
    : config_(RequireElementsOfModel(std::move(config), element_set)),
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
      generator_(MakeGenerator(seed, restart_steps)),
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

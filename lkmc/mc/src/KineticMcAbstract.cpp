#include "KineticMcAbstract.h"
#include "EnergyPredictor.h"
#include <utility>
#include <chrono>
#include <mpi.h>
namespace mc {

KineticMcAbstract::KineticMcAbstract(cfg::Config config,
                                     const unsigned long long int log_dump_steps,
                                     const unsigned long long int config_dump_steps,
                                     const unsigned long long int maximum_steps,
                                     const unsigned long long int thermodynamic_averaging_steps,
                                     const double temperature,
                                     const std::set<Element> &element_set,
                                     const unsigned long long int restart_steps,
                                     const double restart_energy,
                                     const double restart_time,
                                     const std::string &json_coefficients_filename)
    : config_(std::move(config)),
      vacancy_atom_id_(config_.GetVacancyAtomIndex()),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_steps_(maximum_steps),
      beta_(1.0 / constants::kBoltzmann / temperature),
      steps_(restart_steps),
      energy_(restart_energy),
      initial_absolute_energy_(0.0),
      time_(restart_time),
      thermodynamic_averaging_(thermodynamic_averaging_steps),
      energy_predictor_(json_coefficients_filename,
                        config_, element_set, 100000),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      ofs_("lkmc_log.txt", std::ofstream::out | std::ofstream::app) {
  ofs_.precision(16);
  pred::EnergyPredictor total_energy_predictor(json_coefficients_filename, element_set);
  initial_absolute_energy_ = total_energy_predictor.GetEnergy(config_);
}
KineticMcAbstract::~KineticMcAbstract() = default;

void KineticMcAbstract::Dump() const {
  thermodynamic_averaging_.AddEnergy(energy_);
  if (steps_ == 0) {
    config_.WriteLattice("lattice.txt");
    config_.WriteElement("element.txt");
    ofs_ << "steps\ttime\tenergy\taverage_energy\tabsolute_energy\tEa\tdE\ttype" << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteMap("map" + std::to_string(steps_) + ".txt");
  }
  unsigned long long int log_dump_steps;
  if (steps_ > 10 * log_dump_steps_) {
    log_dump_steps = log_dump_steps_;
  } else {
    log_dump_steps = static_cast<unsigned long long int>(
        std::pow(10, static_cast<unsigned long long int>(std::log10(
            steps_ + 1) - 1)));
    log_dump_steps = std::max(log_dump_steps, static_cast<unsigned long long int>(1));
    log_dump_steps = std::min(log_dump_steps, log_dump_steps_);
  }
  if (steps_ % log_dump_steps == 0) {
    ofs_ << steps_ << '\t' << time_ << '\t' << energy_ << '\t'
         << thermodynamic_averaging_.GetThermodynamicAverage(beta_) << '\t'
         << initial_absolute_energy_ + energy_ << '\t'
         << selected_event_.GetForwardBarrier() << '\t' << selected_event_.GetEnergyChange() << '\t'
         << config_.GetElementAtLatticeId(selected_event_.GetIdJumpPair().second).GetString()
         << std::endl;
  }
}
size_t KineticMcAbstract::SelectEvent() const {
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double random_number = distribution(generator_);
  auto it = std::lower_bound(event_k_i_list_.begin(),
                             event_k_i_list_.end(),
                             random_number,
                             [](const auto &lhs, double value) {
                               return lhs.GetCumulativeProvability() < value;
                             });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == event_k_i_list_.cend()) {
    it--;
  }
  return static_cast<size_t>(std::distance(event_k_i_list_.begin(), it));
}
} // mc
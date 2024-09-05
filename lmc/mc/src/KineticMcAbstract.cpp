#include "KineticMcAbstract.h"

namespace mc {

KineticMcFirstAbstract::KineticMcFirstAbstract(cfg::Config config,
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
                                               const Vector_t &vacancy_trajectory,
                                               bool is_early_stop)
    : McAbstract(std::move(config),
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
                 "kmc_log.txt"),
      vacancy_migration_predictor_lru_(json_coefficients_filename, config_, element_set, 10000000),
      time_temperature_interpolator_(time_temperature_filename),
      is_time_temperature_interpolator_(!time_temperature_filename.empty()),
      rate_corrector_(config_.GetVacancyConcentration(), config_.GetSoluteConcentration(Element("Al"))),
      is_rate_corrector_(is_rate_corrector),
      vacancy_lattice_id_(config_.GetVacancyLatticeId()),
      vacancy_trajectory_(vacancy_trajectory),
      is_early_stop_(is_early_stop) {}

KineticMcFirstAbstract::~KineticMcFirstAbstract() = default;

void KineticMcFirstAbstract::UpdateTemperature() {
  if (is_time_temperature_interpolator_) {
    temperature_ = time_temperature_interpolator_.GetTemperature(time_);
    beta_ = 1.0 / constants::kBoltzmann / temperature_;
  }
}

double KineticMcFirstAbstract::GetTimeCorrectionFactor() {
  if (is_rate_corrector_) {
    return rate_corrector_.GetTimeCorrectionFactor(temperature_);
  }
  return 1.0;
}

void KineticMcFirstAbstract::Dump() const {
  if (is_restarted_) {
    is_restarted_ = false;
    return;
  }
  if (world_rank_ != 0) {
    return;
  }
  if (steps_ == 0) {
    config_.WriteLattice("lattice.txt");
    config_.WriteElement("element.txt");
    ofs_ << "steps\ttime\ttemperature\tenergy\tEa\tdE\tselected\tvac1\tvac2\tvac3" << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    // config_.WriteMap("map" + std::to_string(step_) + ".txt");
    config_.WriteConfig(std::to_string(steps_) + ".cfg.gz");
  }
  unsigned long long int log_dump_steps;
  if (steps_ > 10 * log_dump_steps_) {
    log_dump_steps = log_dump_steps_;
  } else {
    log_dump_steps = static_cast<unsigned long long int>(
        std::pow(10, static_cast<unsigned long long int>(std::log10(steps_ + 1) - 1)));
    log_dump_steps = std::max(log_dump_steps, static_cast<unsigned long long int>(1));
    log_dump_steps = std::min(log_dump_steps, log_dump_steps_);
  }
  if (steps_ % log_dump_steps == 0) {
    ofs_ << steps_ << '\t' << time_ << '\t' << temperature_ << '\t' << energy_ << '\t' << event_k_i_.GetForwardBarrier()
         << '\t' << event_k_i_.GetEnergyChange() << '\t'
         << config_.GetAtomIdFromLatticeId(event_k_i_.GetIdJumpPair().second) << '\t' << vacancy_trajectory_
         << std::endl;
  }
}

size_t KineticMcFirstAbstract::SelectEvent() const {
  const double random_number = unit_distribution_(generator_);
  auto it = std::lower_bound(
      event_k_i_list_.begin(), event_k_i_list_.end(), random_number, [](const auto &lhs, double value) {
        return lhs.GetCumulativeProbability() < value;
      });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == event_k_i_list_.cend()) {
    it--;
  }
  if (world_size_ > 1) {
    int event_id = static_cast<int>(it - event_k_i_list_.cbegin());
    MPI_Bcast(&event_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return static_cast<size_t>(event_id);
  } else {
    return static_cast<size_t>(it - event_k_i_list_.cbegin());
  }
}

void KineticMcFirstAbstract::Debug(double one_step_time) const {
  if (world_rank_ == 0) {
    if (std::isnan(one_step_time) or std::isinf(one_step_time) or one_step_time < 0.0) {
      config_.WriteConfig("debug" + std::to_string(steps_) + ".cfg.gz");
      std::cerr << "Invalid time step: " << one_step_time << std::endl;
      std::cerr << "For each event: Energy Barrier, Energy Change, Probability, " << std::endl;
      for (auto &event: event_k_i_list_) {
        std::cerr << event.GetForwardBarrier() << '\t' << event.GetEnergyChange() << '\t'
                  << event.GetCumulativeProbability() << std::endl;
      }
      throw std::runtime_error("Invalid time step");
    }
  }
}

void KineticMcFirstAbstract::OneStepSimulation() {
  UpdateTemperature();
  thermodynamic_averaging_.AddEnergy(energy_);
  BuildEventList();
  double one_step_time = CalculateTime() * GetTimeCorrectionFactor();
  Debug(one_step_time);
  event_k_i_ = event_k_i_list_[SelectEvent()];

  IsEscaped();
  Dump();

  // modify
  time_ += one_step_time;
  energy_ += event_k_i_.GetEnergyChange();
  absolute_energy_ += event_k_i_.GetEnergyChange();

  vacancy_trajectory_ +=
      cfg::GetRelativeDistanceVectorLattice(config_.GetLatticeVector()[vacancy_lattice_id_],
                                            config_.GetLatticeVector()[event_k_i_.GetIdJumpPair().second]) *
      config_.GetBasis();

  config_.LatticeJump(event_k_i_.GetIdJumpPair());
  ++steps_;
  vacancy_lattice_id_ = event_k_i_.GetIdJumpPair().second;
}

void KineticMcFirstAbstract::Simulate() {
  while (steps_ <= maximum_steps_) {
    OneStepSimulation();
  }
}

void KineticMcFirstAbstract::IsEscaped() {
  if (is_early_stop_) {
    size_t Al_count = 0;
    for (auto neighbor_lattice_id: config_.GetFirstNeighborsAdjacencyList()[vacancy_lattice_id_]) {
      if (config_.GetElementAtLatticeId(neighbor_lattice_id) == Element("Al")) {
        Al_count++;
      }
    }
    if (Al_count == constants::kNumFirstNearestNeighbors) {
      if (world_rank_ == 0) {
      config_.WriteConfig("escaped.cfg.gz");
      std::cout << "t_exit: " << time_ << std::endl;
      std::cout << "steps: " << steps_ << std::endl;
      }
      steps_ = maximum_steps_ + 1;
    }
  }
}

KineticMcChainAbstract::KineticMcChainAbstract(cfg::Config config,
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
                                               const Vector_t &vacancy_trajectory,
                                               bool is_early_stop)
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
                             vacancy_trajectory,
                             is_early_stop),
      previous_j_lattice_id_(config_.GetFirstNeighborsAdjacencyList()[vacancy_lattice_id_][0]) {
  MPI_Op_create(DataSum, 1, &mpi_op_);
  DefineStruct(&mpi_datatype_);
}

void KineticMcChainAbstract::OneStepSimulation() {
  KineticMcFirstAbstract::OneStepSimulation();
  previous_j_lattice_id_ = event_k_i_.GetIdJumpPair().first;
}

KineticMcChainAbstract::~KineticMcChainAbstract() {
  MPI_Op_free(&mpi_op_);
  MPI_Type_free(&mpi_datatype_);
}

}    // namespace mc

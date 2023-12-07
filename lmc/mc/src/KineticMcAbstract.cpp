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
                                               const bool is_rate_corrector)
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
      vacancy_migration_predictor_lru_(json_coefficients_filename, config_, element_set, 100000),
      time_temperature_interpolator_(time_temperature_filename),
      is_time_temperature_interpolator_(!time_temperature_filename.empty()),
      rate_corrector_(config_.GetVacancyConcentration(), config_.GetSoluteConcentration(Element("Al"))),
      is_rate_corrector_(is_rate_corrector),
      vacancy_lattice_id_(config_.GetVacancyLatticeId()),
      unwrapped_vacancy_cartesian_coordinate_(config_.GetLatticeVector()[vacancy_lattice_id_].GetCartesianPosition())
{
}
KineticMcFirstAbstract::~KineticMcFirstAbstract() = default;
void KineticMcFirstAbstract::UpdateTemperature()
{
  if (is_time_temperature_interpolator_) {
    temperature_ = time_temperature_interpolator_.GetTemperature(time_);
    beta_ = 1.0 / constants::kBoltzmann / temperature_;
  }
}
double KineticMcFirstAbstract::GetTimeCorrectionFactor()
{
  if (is_rate_corrector_) { return rate_corrector_.GetTimeCorrectionFactor(temperature_); }
  return 1.0;
}
void KineticMcFirstAbstract::Dump() const
{
  if (is_restarted_) {
    is_restarted_ = false;
    return;
  }
  if (world_rank_ != 0) { return; }
  if (steps_ == 0) {
    config_.WriteLattice("lattice.txt");
    config_.WriteElement("element.txt");
    ofs_ << "steps\ttime\ttemperature\tenergy\taverage_energy\tabsolute_energy\tEa\tdE\ttype" << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    // config_.WriteMap("map" + std::to_string(steps_) + ".txt");
    config_.WriteConfig(std::to_string(steps_) + ".cfg.gz");
    ofs_.flush();
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
    ofs_ << steps_ << '\t' << time_ << '\t' << temperature_ << '\t' << energy_ << '\t'
      << thermodynamic_averaging_.GetThermodynamicAverage(beta_) << '\t'
      << absolute_energy_ << '\t'
      << event_k_i_.GetForwardBarrier() << '\t'
      << event_k_i_.GetEnergyChange() << '\t'
      << config_.GetElementAtLatticeId(event_k_i_.GetIdJumpPair().first).GetString()
      << '\n';
    // // unsigned long long step
    // ofs_.write(reinterpret_cast<const char *>(&steps_), sizeof(steps_));
    // // double time
    // ofs_.write(reinterpret_cast<const char *>(&time_), sizeof(time_));
    // // double temperature
    // ofs_.write(reinterpret_cast<const char *>(&temperature_), sizeof(temperature_));
    // // double energy
    // ofs_.write(reinterpret_cast<const char *>(&energy_), sizeof(energy_));
    //
    // // double Ea
    // const double Ea = event_k_i_.GetForwardBarrier();
    // ofs_.write(reinterpret_cast<const char *>(&Ea), sizeof(Ea));
    // // double dE
    // const double dE = event_k_i_.GetEnergyChange();
    // ofs_.write(reinterpret_cast<const char *>(&dE), sizeof(dE));
    //
    // // size_t selected_index
    // const size_t selected_id = config_.GetAtomIdFromLatticeId(event_k_i_.GetIdJumpPair().second);
    // ofs_.write(reinterpret_cast<const char *>(&selected_id), sizeof(selected_id));

    // // type
    // const std::string type = config_.GetElementAtLatticeId(event_k_i_.GetIdJumpPair().second).GetString();
    // const size_t string_size = type.size();
    // ofs_.write(reinterpret_cast<const char *>(&string_size), sizeof(string_size));
    // ofs_.write(type.c_str(), static_cast<int>(string_size));

    // unwraped vacancy position
    // ofs_.write(reinterpret_cast<const char *>(&unwrapped_vacancy_cartesian_coordinate_),
    //            sizeof(unwrapped_vacancy_cartesian_coordinate_));

  }
}
size_t KineticMcFirstAbstract::SelectEvent() const
{
  const double random_number = unit_distribution_(generator_);
  auto it = std::lower_bound(
      event_k_i_list_.begin(), event_k_i_list_.end(), random_number, [](const auto &lhs, double value) {
        return lhs.GetCumulativeProvability() < value;
      });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == event_k_i_list_.cend()) { it--; }
  if (world_size_ > 1) {
    int event_id = static_cast<int>(it - event_k_i_list_.cbegin());
    MPI_Bcast(&event_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return static_cast<size_t>(event_id);
  } else {
    return static_cast<size_t>(it - event_k_i_list_.cbegin());
  }
}
void KineticMcFirstAbstract::OneStepSimulation()
{
  UpdateTemperature();
  thermodynamic_averaging_.AddEnergy(energy_);
  BuildEventList();
  event_k_i_ = event_k_i_list_[SelectEvent()];
  Dump();
  // modify
  time_ += (CalculateTime() * GetTimeCorrectionFactor());
  energy_ += event_k_i_.GetEnergyChange();
  absolute_energy_ += event_k_i_.GetEnergyChange();

  unwrapped_vacancy_cartesian_coordinate_ +=
      cfg::GetRelativeDistanceVectorLattice(config_.GetLatticeVector()[vacancy_lattice_id_],
                                            config_.GetLatticeVector()[event_k_i_.GetIdJumpPair().second]) *
      config_.GetBasis();

  config_.LatticeJump(event_k_i_.GetIdJumpPair());
  ++steps_;
  vacancy_lattice_id_ = event_k_i_.GetIdJumpPair().second;
}
void KineticMcFirstAbstract::Simulate()
{
  while (steps_ <= maximum_steps_) { OneStepSimulation(); }
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
                                               const bool is_rate_corrector)
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
                             is_rate_corrector),
      previous_j_lattice_id_(config_.GetFirstNeighborsAdjacencyList()[vacancy_lattice_id_][0])
{
  MPI_Op_create(DataSum, 1, &mpi_op_);
  DefineStruct(&mpi_datatype_);
}
void KineticMcChainAbstract::OneStepSimulation()
{
  KineticMcFirstAbstract::OneStepSimulation();
  previous_j_lattice_id_ = event_k_i_.GetIdJumpPair().first;
}
KineticMcChainAbstract::~KineticMcChainAbstract()
{
  MPI_Op_free(&mpi_op_);
  MPI_Type_free(&mpi_datatype_);
}

}    // namespace mc

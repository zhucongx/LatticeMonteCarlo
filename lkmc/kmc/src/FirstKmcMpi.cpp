#include "FirstKmcMpi.h"
#include <utility>
#include <chrono>
namespace kmc {
FirstKmcMpi::FirstKmcMpi(cfg::Config config,
                         unsigned long long int log_dump_steps,
                         unsigned long long int config_dump_steps,
                         unsigned long long int maximum_number,
                         double temperature,
                         const std::set<Element> &type_set,
                         unsigned long long int restart_steps,
                         double restart_energy,
                         double restart_time,
                         const std::string &json_coefficients_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      beta_(1 / kBoltzmannConstant / temperature),
      steps_(restart_steps),
      energy_(restart_energy),
      time_(restart_time),
      vacancy_index_(cfg::GetVacancyAtomIndex(config_)),
      energy_predictor_(json_coefficients_filename,
                        config_, type_set, 100000),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
  event_list_.resize(kEventListSize);
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_size != kEventListSize) {
    std::cout << "Must use " << kEventListSize << " precesses. Terminating.\n";
    MPI_Finalize();
    exit(0);
  }
  if (mpi_rank_ == 0) {
    std::cout << "Using " << mpi_size << " processes." << std::endl;
  }
}
FirstKmcMpi::~FirstKmcMpi() {
  MPI_Finalize();
}
void FirstKmcMpi::Dump(std::ofstream &ofs) const {
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << time_ << '\t' << energy_ << '\t' << one_step_barrier_ << '\t'
        << one_step_energy_change_ << '\t' << migrating_element_.GetString() << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteCfg(std::to_string(steps_) + ".cfg", true);
  }
}
void FirstKmcMpi::BuildEventListParallel() {
  total_rate_ = 0;
  const auto neighbor_index =
      config_.GetFirstNeighborsAtomIdVectorOfAtom(vacancy_index_)[static_cast<size_t>(mpi_rank_)];
  JumpEvent event
      ({vacancy_index_, neighbor_index},
       energy_predictor_.GetBarrierAndDiffFromAtomIdPair(config_,
                                                         {vacancy_index_, neighbor_index}),
       beta_);
  const double this_rate = event.GetForwardRate();

  MPI_Allgather(&event,
                sizeof(JumpEvent),
                MPI_BYTE,
                event_list_.data(),
                sizeof(JumpEvent),
                MPI_BYTE,
                MPI_COMM_WORLD);
  MPI_Allreduce(&this_rate, &total_rate_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double cumulative_probability = 0.0;
  for (auto &event_it: event_list_) {
    event_it.CalculateProbability(total_rate_);
    cumulative_probability += event_it.GetProbability();
    event_it.SetCumulativeProbability(cumulative_probability);
  }
}
size_t FirstKmcMpi::SelectEvent() const {
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);

  const double random_number = distribution(generator_);
  auto it = std::lower_bound(event_list_.begin(),
                             event_list_.end(),
                             random_number,
                             [](const auto &lhs, double value) {
                               return lhs.GetCumulativeProvability() < value;
                             });
  // If not find (maybe generated 1), which rarely happens, returns the last event
  if (it == event_list_.cend()) {
    it--;
  }
  return static_cast<size_t>(std::distance(event_list_.begin(), it));
}
void FirstKmcMpi::Simulate() {
  std::ofstream ofs("kmc_log.txt", std::ofstream::out | std::ofstream::app);
  if (mpi_rank_ == 0) {
    ofs << "steps\ttime\tenergy\tEa\tdE\ttype\n";
    ofs.precision(8);
  }

  while (steps_ < maximum_number_) {
    if (mpi_rank_ == 0) {
      Dump(ofs);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    BuildEventListParallel();
    JumpEvent selected_event;
    if (mpi_rank_ == 0) {
      selected_event = event_list_[SelectEvent()];
    }
    MPI_Bcast(&selected_event, sizeof(JumpEvent), MPI_BYTE, 0, MPI_COMM_WORLD);

    atom_id_jump_pair_ = selected_event.GetAtomIdJumpPair();

    // if (!CheckAndSolveEquilibrium(ofs)) {
    // update time and energy
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    one_step_time_change_ = -log(distribution(generator_)) / total_rate_ / kPrefactor;
    time_ += one_step_time_change_;
    one_step_energy_change_ = selected_event.GetEnergyChange();
    energy_ += one_step_energy_change_;
    one_step_barrier_ = selected_event.GetForwardBarrier();
    migrating_element_ = config_.GetAtomVector().at(atom_id_jump_pair_.second).GetElement();
    config_.AtomJump(atom_id_jump_pair_);
    // }
    ++steps_;
  }
}
} // namespace kmc
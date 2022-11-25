#include "KineticMcChainOmpi.h"
#include "EnergyPredictor.h"
namespace mc {
//  j -> k -> i ->l
//       |
// current position
struct MpiData {
  double beta_bar_k{0.0};
  double beta_k{0.0};
  double gamma_bar_k_j{0.0};
  double gamma_k_j{0.0};
  double beta_k_j{0.0};
  double alpha_k_j{0.0};
  double ts_numerator{0.0};
  double ts_j_numerator{0.0};
};

void DataSum(void *input_buffer, void *output_buffer, int *len, MPI_Datatype *datatype) {
  auto *input = static_cast<MpiData *>(input_buffer);
  auto *output = static_cast<MpiData *>(output_buffer);
  for (int i = 0; i < *len; ++i) {
    output[i].beta_bar_k += input[i].beta_bar_k;
    output[i].beta_k += input[i].beta_k;
    output[i].gamma_bar_k_j += input[i].gamma_bar_k_j;
    output[i].gamma_k_j += input[i].gamma_k_j;
    output[i].beta_k_j += input[i].beta_k_j;
    output[i].alpha_k_j += input[i].alpha_k_j;
    output[i].ts_numerator += input[i].ts_numerator;
    output[i].ts_j_numerator += input[i].ts_j_numerator;
  }
}
void DefineStruct(MPI_Datatype *datatype) {
  const int count = 8;
  int block_lens[count];
  MPI_Datatype types[count];
  MPI_Aint displacements[count];

  for (int i = 0; i < count; i++) {
    types[i] = MPI_DOUBLE;
    block_lens[i] = 1;
  }
  displacements[0] = offsetof(MpiData, beta_bar_k);
  displacements[1] = offsetof(MpiData, beta_k);
  displacements[2] = offsetof(MpiData, gamma_bar_k_j);
  displacements[3] = offsetof(MpiData, gamma_k_j);
  displacements[4] = offsetof(MpiData, beta_k_j);
  displacements[5] = offsetof(MpiData, alpha_k_j);
  displacements[6] = offsetof(MpiData, ts_numerator);
  displacements[7] = offsetof(MpiData, ts_j_numerator);

  MPI_Type_create_struct(count, block_lens, displacements, types, datatype);
  MPI_Type_commit(datatype);
}

KineticMcChainOmpi::KineticMcChainOmpi(cfg::Config config,
                                       unsigned long long int log_dump_steps,
                                       unsigned long long int config_dump_steps,
                                       unsigned long long int maximum_number,
                                       double temperature,
                                       const std::set<Element> &element_set,
                                       unsigned long long int restart_steps,
                                       double restart_energy,
                                       double restart_time,
                                       const std::string &json_coefficients_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_number_(maximum_number),
      beta_(1.0 / constants::kBoltzmann / temperature),
      steps_(restart_steps),
      energy_(restart_energy),
      time_(restart_time),
      vacancy_index_(config_.GetVacancyAtomIndex()),
      previous_j_(config_.GetFirstNeighborsAtomIdVectorOfAtom(vacancy_index_)[0]),
      energy_predictor_(json_coefficients_filename,
                        config_, element_set, 100000),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())) {
  event_k_i_list_.resize(kEventListSize);
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_size != kEventListSize) {
    std::cout << "Must use " << kEventListSize << " precesses. Terminating...\n" << std::endl;
    MPI_Finalize();
    exit(0);
  }
  if (world_rank_ == 0) {
    std::cout << "Using " << mpi_size << " processes." << std::endl;
#pragma omp parallel master default(none) shared(std::cout)
    {
      std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
    }
  }
  MPI_Op_create(DataSum, 1, &mpi_op_);
  DefineStruct(&mpi_datatype_);
  if (world_rank_ == 0) {
    std::ofstream ofs("lkmc_log.txt", std::ofstream::out | std::ofstream::app);
    ofs.precision(16);

    pred::EnergyPredictor total_energy_predictor(json_coefficients_filename, element_set);
    ofs << "initial_energy = " << total_energy_predictor.GetEnergy(config_) << std::endl;
    ofs << "steps\ttime\tenergy\tEa\tdE\ttype" << std::endl;
  }
}
KineticMcChainOmpi::~KineticMcChainOmpi() {
  MPI_Op_free(&mpi_op_);
  MPI_Finalize();
}
void KineticMcChainOmpi::Dump(std::ofstream &ofs) const {
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << time_ << '\t' << energy_ << '\t' << one_step_barrier_ << '\t'
        << one_step_energy_change_ << '\t' << migrating_element_.GetString() << std::endl;
  }
  if (steps_ == 0) {
    config_.WriteLattice("lattice.txt");
    config_.WriteElement("element.txt");
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteMap("map" + std::to_string(steps_) + ".txt");
  }
}
// update  first_event_ki and l_index_list for each process
void KineticMcChainOmpi::BuildFirstEventKIAndGetTotalRates() {
  const auto i_index = config_.GetFirstNeighborsAtomIdVectorOfAtom(
      vacancy_index_)[static_cast<size_t>(world_rank_)];
  size_t it = 0;
  for (auto l_index: config_.GetFirstNeighborsAtomIdVectorOfAtom(i_index)) {
    if (l_index == vacancy_index_) {
      l_index = i_index;
    }
    l_index_list_[it] = l_index;
    ++it;
  }

  total_rate_k_ = 0.0;
  total_rate_i_ = 0.0;

  config_.AtomJump({vacancy_index_, i_index});
#pragma omp parallel for default(none) shared(i_index) reduction(+:total_rate_i_)
  for (size_t ii = 0; ii < kEventListSize; ++ii) {
    const auto l_index = l_index_list_[ii];
    JumpEvent event_i_l({vacancy_index_, l_index},
                        energy_predictor_.GetBarrierAndDiffFromAtomIdPair(
                            config_, {vacancy_index_, l_index}),
                        beta_);
    if (l_index == i_index) {
      event_k_i_ = event_i_l.GetReverseJumpEvent();
    }
    auto r_i_l = event_i_l.GetForwardRate();
    total_rate_i_ += r_i_l;
  }
  config_.AtomJump({vacancy_index_, i_index});

  double rate_i_k = event_k_i_.GetForwardRate();
  MPI_Allreduce(&rate_i_k, &total_rate_k_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  event_k_i_.CalculateProbability(total_rate_k_);
}

double KineticMcChainOmpi::UpdateIndirectProbabilityAndCalculateTime() {
  const auto probability_k_i = event_k_i_.GetProbability();
  const auto probability_i_k = event_k_i_.GetBackwardRate() / total_rate_i_;
  const double beta_bar_k_i = probability_k_i * probability_i_k;
  const double beta_k_i = probability_k_i * (1 - probability_i_k);
  bool is_previous_event = event_k_i_.GetAtomIdJumpPair().second == previous_j_;

  // Time in first order Kmc, same for all
  const double t_1 = 1 / total_rate_k_ / constants::kPrefactor;
  // Time in neighbors, different for each process
  const double t_i = 1 / total_rate_i_ / constants::kPrefactor;
  MpiData mpi_data_helper{
      beta_bar_k_i,
      beta_k_i,
      is_previous_event ? 0.0 : beta_bar_k_i,
      is_previous_event ? 0.0 : beta_k_i,
      is_previous_event ? beta_k_i : 0.0,
      is_previous_event ? probability_k_i : 0.0,
      (t_1 + t_i) * beta_bar_k_i,
      is_previous_event ? 0.0 : (t_1 + t_i) * beta_bar_k_i
  };
  MpiData mpi_data{};
  MPI_Allreduce(&mpi_data_helper, &mpi_data, 1, mpi_datatype_, mpi_op_, MPI_COMM_WORLD);

  const auto beta_bar_k = mpi_data.beta_bar_k;
  const auto beta_k = mpi_data.beta_k;
  const auto gamma_bar_k_j = mpi_data.gamma_bar_k_j;
  const auto gamma_k_j = mpi_data.gamma_k_j;
  const auto beta_k_j = mpi_data.beta_k_j;
  const auto alpha_k_j = mpi_data.alpha_k_j;
  const auto ts = mpi_data.ts_numerator / beta_bar_k;
  const auto ts_j = mpi_data.ts_j_numerator / gamma_bar_k_j;

  const double one_over_one_minus_a_j = 1 / (1 - alpha_k_j);
  const double t_2 = one_over_one_minus_a_j
      * (gamma_k_j * t_1 + gamma_bar_k_j * (ts_j + t_1 + beta_bar_k / beta_k * ts));

  double indirect_probability;
  if (is_previous_event) {
    indirect_probability = one_over_one_minus_a_j * (gamma_bar_k_j / beta_k) * beta_k_j;
  } else {
    indirect_probability = one_over_one_minus_a_j * (1 + gamma_bar_k_j / beta_k) * beta_k_i;
  }
  event_k_i_.SetProbability(indirect_probability);
  MPI_Allgather(&event_k_i_,
                sizeof(JumpEvent),
                MPI_BYTE,
                event_k_i_list_.data(),
                sizeof(JumpEvent),
                MPI_BYTE,
                MPI_COMM_WORLD);
  double cumulative_probability = 0.0;
  for (auto &event: event_k_i_list_) {
    cumulative_probability += event.GetProbability();
    event.SetCumulativeProbability(cumulative_probability);
  }
  return t_2;
}

// run this on this first process
size_t KineticMcChainOmpi::SelectEvent() const {
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

void KineticMcChainOmpi::Simulate() {
  std::ofstream ofs("lkmc_log.txt", std::ofstream::out | std::ofstream::app);
  if (world_rank_ == 0) {
    ofs.precision(16);
  }
  while (steps_ <= maximum_number_) {
    if (world_rank_ == 0) {
      Dump(ofs);
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    BuildFirstEventKIAndGetTotalRates();
    one_step_time_change_ = UpdateIndirectProbabilityAndCalculateTime();
    time_ += one_step_time_change_;

    JumpEvent selected_event;
    if (world_rank_ == 0) {
      selected_event = event_k_i_list_[SelectEvent()];
    }
    MPI_Bcast(&selected_event, sizeof(JumpEvent), MPI_BYTE, 0, MPI_COMM_WORLD);
    atom_id_jump_pair_ = selected_event.GetAtomIdJumpPair();

    one_step_energy_change_ = selected_event.GetEnergyChange();
    energy_ += one_step_energy_change_;
    one_step_barrier_ = selected_event.GetForwardBarrier();
    migrating_element_ = config_.GetAtomVector().at(atom_id_jump_pair_.second).GetElement();
    config_.AtomJump(atom_id_jump_pair_);
    previous_j_ = atom_id_jump_pair_.second;
    ++steps_;
  }
}

} // mc
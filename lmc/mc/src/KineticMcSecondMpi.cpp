#include "KineticMcSecondMpi.h"
namespace mc {
namespace {
void DataSum(void *input_buffer,
             void *output_buffer,
             int *len,
             [[maybe_unused]] MPI_Datatype *datatype) {
  auto *input = static_cast<SecondOrderSums *>(input_buffer);
  auto *output = static_cast<SecondOrderSums *>(output_buffer);
  for (int i = 0; i < *len; ++i) {
    output[i] += input[i];
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
  displacements[0] = offsetof(SecondOrderSums, beta_bar_k);
  displacements[1] = offsetof(SecondOrderSums, beta_k);
  displacements[2] = offsetof(SecondOrderSums, gamma_bar_k_j);
  displacements[3] = offsetof(SecondOrderSums, gamma_k_j);
  displacements[4] = offsetof(SecondOrderSums, beta_k_j);
  displacements[5] = offsetof(SecondOrderSums, alpha_k_j);
  displacements[6] = offsetof(SecondOrderSums, ts_numerator);
  displacements[7] = offsetof(SecondOrderSums, ts_j_numerator);
  MPI_Type_create_struct(count, block_lens, displacements, types, datatype);
  MPI_Type_commit(datatype);
}
}    // namespace
//  j -> k -> i ->l
//       |
// current position
KineticMcSecondMpi::KineticMcSecondMpi(cfg::Config config,
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
  if (world_size_ != kEventListSize) {
    if (world_rank_ == 0) {
      std::cout << "Must use " << kEventListSize << " processes. Terminating...\n" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (world_rank_ == 0) {
    std::cout << "Using " << world_size_ << " processes." << std::endl;
  }
  MPI_Op_create(DataSum, 1, &mpi_op_);
  DefineStruct(&mpi_datatype_);
}

KineticMcSecondMpi::~KineticMcSecondMpi() {
  MPI_Op_free(&mpi_op_);
  MPI_Type_free(&mpi_datatype_);
}

// update first_event_ki for each process
void KineticMcSecondMpi::BuildEventList() {
  const auto k_lattice_id = vacancy_lattice_id_;
  const auto i_lattice_id =
      config_.GetFirstNeighborsAdjacencyList()[k_lattice_id][static_cast<size_t>(world_rank_)];

  total_rate_k_ = 0.0;

  config_.LatticeJump({k_lattice_id, i_lattice_id});
  EvaluateSecondJumps(k_lattice_id, i_lattice_id);
  config_.LatticeJump({i_lattice_id, k_lattice_id});

  double rate_k_i = event_k_i_.GetForwardRate();
  MPI_Allreduce(&rate_k_i, &total_rate_k_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  event_k_i_.CalculateProbability(total_rate_k_);
}

void KineticMcSecondMpi::EvaluateSecondJumps(const size_t k_lattice_id, const size_t i_lattice_id) {
  total_rate_i_ = 0.0;
  for (const auto l_lattice_id: config_.GetFirstNeighborsAdjacencyList()[i_lattice_id]) {
    JumpEvent event_i_l({i_lattice_id, l_lattice_id},
                        vacancy_migration_predictor_lru_.GetBarrierAndDiffFromLatticeIdPair(
                            config_, {i_lattice_id, l_lattice_id}),
                        beta_);
    if (l_lattice_id == k_lattice_id) {
      event_k_i_ = event_i_l.GetReverseJumpEvent();
    }
    total_rate_i_ += event_i_l.GetForwardRate();
  }
}

double KineticMcSecondMpi::CalculateTime() {
  const auto summands = GetSecondOrderSummands(event_k_i_, total_rate_i_);
  SecondOrderSums sums{};
  MPI_Allreduce(&summands, &sums, 1, mpi_datatype_, mpi_op_, MPI_COMM_WORLD);

  event_k_i_.SetProbability(GetSecondOrderProbability(event_k_i_, summands, sums));
  MPI_Allgather(&event_k_i_,
                sizeof(JumpEvent),
                MPI_BYTE,
                event_k_i_list_.data(),
                sizeof(JumpEvent),
                MPI_BYTE,
                MPI_COMM_WORLD);
  // TODO(perf): Reduce MPI payload by gathering only minimal fields (probability,
  // forward/backward rates or the chosen index) and reconstructing locally.
  double cumulative_probability = 0.0;
  for (auto &event : event_k_i_list_) {
    cumulative_probability += event.GetProbability();
    event.SetCumulativeProbability(cumulative_probability);
  }
  return GetSecondOrderTime(sums);
}
} // mc

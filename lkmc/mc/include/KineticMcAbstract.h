#ifndef LKMC_LKMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#define LKMC_LKMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "JumpEvent.h"
#include "ThermodynamicAveraging.h"
#include "VacancyMigrationPredictorQuarticLru.h"
namespace mc {

class KineticMcAbstract {
    // it is lattice id jump pair based
  public:
    KineticMcAbstract(cfg::Config config,
                      unsigned long long int log_dump_steps,
                      unsigned long long int config_dump_steps,
                      unsigned long long int maximum_steps,
                      unsigned long long int thermodynamic_averaging_steps,
                      double temperature,
                      const std::set<Element> &element_set,
                      unsigned long long int restart_steps,
                      double restart_energy,
                      double restart_time,
                      const std::string &json_coefficients_filename);
    virtual ~KineticMcAbstract();
    virtual void Simulate() = 0;
  protected:
    void Dump() const;
    size_t SelectEvent() const;

    // constants
    static constexpr size_t kEventListSize = constants::kNumFirstNearestNeighbors;

    // config
    cfg::Config config_;
    const size_t vacancy_atom_id_;
    // simulation parameters
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_steps_;
    const double beta_;
    // simulation statistics
    unsigned long long int steps_;
    double energy_;
    double initial_absolute_energy_;
    double time_;
    // helpful properties
    mutable ThermodynamicAveraging thermodynamic_averaging_;
    const pred::VacancyMigrationPredictorQuarticLru energy_predictor_;
    mutable std::mt19937_64 generator_;
    mutable std::ofstream ofs_;

    // simulation variables
    std::array<JumpEvent, kEventListSize> event_k_i_list_{};
    JumpEvent selected_event_k_i_{};
    double total_rate_k_{0.0}; // k would be same for all
};
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

inline void DataSum(void *input_buffer,
                    void *output_buffer,
                    int *len,
                    [[maybe_unused]] MPI_Datatype *datatype) {
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
inline void DefineStruct(MPI_Datatype *datatype) {
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

} // mc

#endif //LKMC_LKMC_MC_INCLUDE_KINETICMCABSTRACT_H_

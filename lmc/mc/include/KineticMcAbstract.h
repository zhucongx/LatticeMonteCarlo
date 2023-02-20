#ifndef LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "McAbstract.h"
#include "JumpEvent.h"
#include "VacancyMigrationPredictorQuarticLru.h"
#include "TimeTemperatureInterpolator.h"
#include "RateCorrector.hpp"
namespace mc {

class KineticMcFirstAbstract : public McAbstract {
    // it is lattice id jump pair based
  public:
    KineticMcFirstAbstract(cfg::Config config,
                           unsigned long long int log_dump_steps,
                           unsigned long long int config_dump_steps,
                           unsigned long long int maximum_steps,
                           unsigned long long int thermodynamic_averaging_steps,
                           unsigned long long int restart_steps,
                           double restart_energy,
                           double restart_time,
                           double temperature,
                           const std::set<Element> &element_set,
                           const std::string &json_coefficients_filename,
                           const std::string &time_temperature_filename,
                           bool is_rate_corrector);
    ~KineticMcFirstAbstract() override;
    KineticMcFirstAbstract(const KineticMcFirstAbstract &) = delete;
    void operator=(const mc::KineticMcFirstAbstract &) = delete;
    void Simulate() override;
  protected:
    void UpdateTemperature();
    double GetTimeCorrectionFactor();
    virtual void Dump() const;
    size_t SelectEvent() const;
    virtual void BuildEventList() = 0;
    virtual double CalculateTime() = 0;
    virtual void OneStepSimulation();
    // constants
    static constexpr size_t kEventListSize = constants::kNumFirstNearestNeighbors;

    // helpful properties
    const pred::VacancyMigrationPredictorQuarticLru vacancy_migration_predictor_lru_;
    const pred::TimeTemperatureInterpolator time_temperature_interpolator_;
    const bool is_time_temperature_interpolator_;
    const pred::RateCorrector rate_corrector_;
    const bool is_rate_corrector_;
    size_t vacancy_lattice_id_;
    std::array<JumpEvent, kEventListSize> event_k_i_list_{};
    JumpEvent event_k_i_{};
    double total_rate_k_{0.0}; // k would be same for all
};
class KineticMcChainAbstract : public KineticMcFirstAbstract {
  public:
    KineticMcChainAbstract(cfg::Config config,
                           unsigned long long int log_dump_steps,
                           unsigned long long int config_dump_steps,
                           unsigned long long int maximum_steps,
                           unsigned long long int thermodynamic_averaging_steps,
                           unsigned long long int restart_steps,
                           double restart_energy,
                           double restart_time,
                           double temperature,
                           const std::set<Element> &element_set,
                           const std::string &json_coefficients_filename,
                           const std::string &time_temperature_filename,
                           bool is_rate_corrector);
    ~KineticMcChainAbstract() override;
    KineticMcChainAbstract(const KineticMcChainAbstract &) = delete;
    void operator=(const mc::KineticMcChainAbstract &) = delete;
  protected:
    void OneStepSimulation() override;
    // helpful properties
    size_t previous_j_lattice_id_;
    double total_rate_i_{0.0}; // i would be different
    std::array<size_t, kEventListSize> l_lattice_id_list_{};

    MPI_Op mpi_op_{};
    MPI_Datatype mpi_datatype_{};
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

#endif //LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_

#ifndef LMC_LMC_MC_INCLUDE_KINETICMCSECONDMPI_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCSECONDMPI_H_
#include <mpi.h>
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
namespace mc {
//  j -> k -> i -> l
//       |
// current position

// Second-order KMC on twelve processes, one for each first jump k -> i. A process evaluates the twelve jumps i -> l that
// can follow its first jump one after the other.
class KineticMcSecondMpi : public KineticMcSecondAbstract {
  public:
    KineticMcSecondMpi(cfg::Config config,
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
                       bool is_early_stop,
                       bool is_solute_disp,
                       unsigned long long int seed);
    ~KineticMcSecondMpi() override;
  protected:
    void BuildEventList() override;
    double CalculateTime() override;
    // Evaluates the twelve jumps i -> l of the vacancy, which the caller has moved to site i. Sets event_k_i_ to the
    // reverse of the jump back to k, and total_rate_i_ to the sum of the twelve rates.
    virtual void EvaluateSecondJumps(size_t k_lattice_id, size_t i_lattice_id);
    // The sum of the rates of the twelve jumps i -> l that can follow the first jump of this process.
    double total_rate_i_{0.0};
    // Add up the summands of the twelve processes.
    MPI_Op mpi_op_{};
    MPI_Datatype mpi_datatype_{};
};

} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCSECONDMPI_H_

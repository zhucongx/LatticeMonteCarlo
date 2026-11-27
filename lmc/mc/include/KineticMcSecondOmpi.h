#ifndef LMC_LMC_MC_INCLUDE_KINETICMCSECONDOMPI_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCSECONDOMPI_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "VacancyMigrationPredictorQuarticLru.h"
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
namespace mc {
//  j -> k -> i -> l
//       |
// current position

class KineticMcSecondOmpi : public KineticMcSecondAbstract {
  public:
    KineticMcSecondOmpi(cfg::Config config,
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
    ~KineticMcSecondOmpi() override;
  protected:
    void BuildEventList() override;
    double CalculateTime() override;
    // The sum of the rates of the twelve jumps i -> l that can follow the first jump of this process.
    double total_rate_i_{0.0};
    // Add up the summands of the twelve processes.
    MPI_Op mpi_op_{};
    MPI_Datatype mpi_datatype_{};
};

} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCSECONDOMPI_H_

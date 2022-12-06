#ifndef LMC_LMC_MC_INCLUDE_MCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_MCABSTRACT_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "Config.h"
#include "ThermodynamicAveraging.h"
namespace mc {

class McAbstract {
  public:
    McAbstract(cfg::Config config,
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
               const std::string &log_filename);
    ~McAbstract();
    virtual void Simulate() = 0;
  protected:
    // config
    cfg::Config config_;
    // simulation parameters
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_steps_;
    // simulation statistics
    unsigned long long int steps_;
    double energy_;
    double initial_absolute_energy_;
    double time_;
    double temperature_;
    double beta_;
    bool is_restarted_;
    // helpful properties
    ThermodynamicAveraging thermodynamic_averaging_;
    mutable std::mt19937_64 generator_;
    mutable std::uniform_real_distribution<double> unit_distribution_;
    mutable std::ofstream ofs_;

    int world_rank_{-1};
    int world_size_{-1};
  protected:

};

} // mc

#endif //LMC_LMC_MC_INCLUDE_MCABSTRACT_H_

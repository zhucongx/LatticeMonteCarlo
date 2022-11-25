#ifndef LKMC_LKMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#define LKMC_LKMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#include <random>
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
    JumpEvent selected_event_{};
    double total_rate_k_{0.0}; // k would be same for all
};

} // mc

#endif //LKMC_LKMC_MC_INCLUDE_KINETICMCABSTRACT_H_

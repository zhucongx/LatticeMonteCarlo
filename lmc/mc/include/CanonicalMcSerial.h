#ifndef LMC_LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
#define LMC_LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "CanonicalMcAbstract.h"
#include "EnergyChangePredictorPairAll.h"
namespace mc {
class CanonicalMcSerial : public CanonicalMcAbstract {
  public:
    CanonicalMcSerial(cfg::Config config,
                      unsigned long long int log_dump_steps,
                      unsigned long long int config_dump_steps,
                      unsigned long long int maximum_steps,
                      unsigned long long int thermodynamic_averaging_steps,
                      double initial_temperature,
                      double decrement_temperature,
                      const std::set<Element> &element_set,
                      const std::string &json_coefficients_filename);
    void Simulate() override;
};
} // mc

#endif //LMC_LMC_MC_INCLUDE_CANONICALMCSERIAL_H_

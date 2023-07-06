/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 4:05 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_MC_INCLUDE_CANONICALMCOMP_H_
#define LMC_MC_INCLUDE_CANONICALMCOMP_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "CanonicalMcAbstract.h"
#include "EnergyChangePredictorPairAll.h"
namespace mc {

class CanonicalMcOmp : public CanonicalMcAbstract {
 public:
  CanonicalMcOmp(cfg::Config config,
                 unsigned long long int log_dump_steps,
                 unsigned long long int config_dump_steps,
                 unsigned long long int maximum_steps,
                 unsigned long long int thermodynamic_averaging_steps,
                 unsigned long long int restart_steps,
                 double restart_energy,
                 double temperature,
      // double initial_temperature,
      // double decrement_temperature,
                 const std::set<Element> &element_set,
                 const std::string &json_coefficients_filename);
  void Simulate() override;
 private:
  void BuildEventVector();
  // helpful properties
  std::vector<std::pair<std::pair<size_t, size_t>, double> > event_vector_{};
  std::unordered_set<size_t> unavailable_position_{};
  size_t num_threads_{};
};

} // mc

#endif //LMC_MC_INCLUDE_CANONICALMCOMP_H_

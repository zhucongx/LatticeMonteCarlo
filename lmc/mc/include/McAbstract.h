/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/26/23 9:26 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_MC_INCLUDE_MCABSTRACT_H_
#define LMC_MC_INCLUDE_MCABSTRACT_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "Config.h"
#include "ThermodynamicAveraging.h"
class McAbstract {
 public:
  McAbstract(Config config,
             unsigned long long int log_dump_steps,
             unsigned long long int config_dump_steps,
             unsigned long long int maximum_steps,
      // unsigned long long int thermodynamic_averaging_steps,
             unsigned long long int restart_steps,
             double restart_energy,
             double restart_time,
             double temperature,
             const std::set<Element> &element_set,
             const std::string &json_coefficients_filename,
             const std::string &log_filename);
  virtual ~McAbstract();
  McAbstract(const McAbstract &) = delete;
  void operator=(const McAbstract &) = delete;
  virtual void Simulate() = 0;
 protected:
  // config
  Config config_;
  // simulation parameters
  const unsigned long long int log_dump_steps_;
  const unsigned long long int config_dump_steps_;
  const unsigned long long int maximum_steps_;
  // simulation statistics
  unsigned long long int steps_;
  double energy_;
  double absolute_energy_;
  double time_;
  double temperature_;
  double beta_;
  mutable bool is_restarted_;
  // helpful properties
  // ThermodynamicAveraging thermodynamic_averaging_;
  mutable std::mt19937_64 generator_;
  mutable std::uniform_real_distribution<double> unit_distribution_;
  mutable std::ofstream ofs_;

  int world_rank_{-1};
  int world_size_{-1};
 protected:

};

#endif //LMC_MC_INCLUDE_MCABSTRACT_H_

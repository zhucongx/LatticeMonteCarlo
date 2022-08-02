#ifndef LKMC_LKMC_API_INCLUDE_HOME_H_
#define LKMC_LKMC_API_INCLUDE_HOME_H_
#include "Parameter.h"
#include "FirstKmcMpi.h"
#include "ChainKmcMpi.h"
#include "Iterator.h"
#include "SimulatedAnnealing.h"
namespace api {
kmc::FirstKmcMpi BuildFirstKmcMpiFromParameter(const Parameter &parameter);
kmc::ChainKmcMpi BuildChainKmcMpiFromParameter(const Parameter &parameter);
ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter);
ansys::Iterator BuildIteratorFromParameter(const Parameter &parameter);
} // namespace api

#endif //LKMC_LKMC_API_INCLUDE_HOME_H_

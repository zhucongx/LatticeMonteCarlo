#ifndef LKMC_LKMC_API_INCLUDE_HOME_H_
#define LKMC_LKMC_API_INCLUDE_HOME_H_
#include "Parameter.h"
#include "FirstKmcMpi.h"
#include "ChainKmcMpi.h"
#include "ChainKmcOmp.h"
#include "ChainKmcOmpi.h"
#include "SimulatedAnnealing.h"
#include "CanonicalMc.h"
#include "CanonicalMcStepT.h"
#include "Iterator.h"

namespace api {
kmc::FirstKmcMpi BuildFirstKmcMpiFromParameter(const Parameter &parameter);
kmc::ChainKmcMpi BuildChainKmcMpiFromParameter(const Parameter &parameter);
kmc::ChainKmcOmp BuildChainKmcOmpFromParameter(const Parameter &parameter);
kmc::ChainKmcOmpi BuildChainKmcOmpiFromParameter(const Parameter &parameter);
ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter);
ansys::CanonicalMc BuildCanonicalMcFromParameter(const Parameter &parameter);
ansys::CanonicalMcStepT BuildCanonicalMcStepTFromParameter(const Parameter &parameter);
ansys::Iterator BuildIteratorFromParameter(const Parameter &parameter);
} // namespace api

#endif //LKMC_LKMC_API_INCLUDE_HOME_H_

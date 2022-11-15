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
mc::FirstKmcMpi BuildFirstKmcMpiFromParameter(const Parameter &parameter);
mc::ChainKmcMpi BuildChainKmcMpiFromParameter(const Parameter &parameter);
mc::ChainKmcOmp BuildChainKmcOmpFromParameter(const Parameter &parameter);
mc::ChainKmcOmpi BuildChainKmcOmpiFromParameter(const Parameter &parameter);
ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter);
mc::CanonicalMc BuildCanonicalMcFromParameter(const Parameter &parameter);
mc::CanonicalMcStepT BuildCanonicalMcStepTFromParameter(const Parameter &parameter);
ansys::Iterator BuildIteratorFromParameter(const Parameter &parameter);
} // namespace api

#endif //LKMC_LKMC_API_INCLUDE_HOME_H_

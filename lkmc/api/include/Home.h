#ifndef LKMC_LKMC_API_INCLUDE_HOME_H_
#define LKMC_LKMC_API_INCLUDE_HOME_H_
#include "Parameter.h"
#include "KineticMcFirstMpi.h"
#include "KineticMcChainMpi.h"
#include "KineticMcChainOmp.h"
#include "KineticMcChainOmpi.h"
#include "SimulatedAnnealing.h"
#include "CanonicalMc.h"
#include "CanonicalMcStepT.h"
#include "Iterator.h"

namespace api {
mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter(const Parameter &parameter);
mc::KineticMcChainMpi BuildKineticMcChainMpiFromParameter(const Parameter &parameter);
mc::KineticMcChainOmp BuildKineticMcChainOmpFromParameter(const Parameter &parameter);
mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter &parameter);
ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter);
mc::CanonicalMc BuildCanonicalMcFromParameter(const Parameter &parameter);
mc::CanonicalMcStepT BuildCanonicalMcStepTFromParameter(const Parameter &parameter);
ansys::Iterator BuildIteratorFromParameter(const Parameter &parameter);
} // namespace api

#endif //LKMC_LKMC_API_INCLUDE_HOME_H_

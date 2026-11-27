#ifndef LMC_LMC_API_INCLUDE_HOME_H_
#define LMC_LMC_API_INCLUDE_HOME_H_
#include "Parameter.h"
#include "KineticMcFirstMpi.h"
#include "KineticMcFirstSerial.h"
#include "KineticMcSecondMpi.h"
#include "KineticMcSecondOmpi.h"
#include "SimulatedAnnealing.h"
#include "CanonicalMcSerial.h"
#include "CanonicalMcOmp.h"
#include "Traverse.h"

namespace api {
void Print(const Parameter &parameter);
void Run(const Parameter &parameter);
mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter(const Parameter &parameter);
mc::KineticMcFirstSerial BuildKineticMcFirstSerialFromParameter(const Parameter &parameter);
mc::KineticMcSecondMpi BuildKineticMcSecondMpiFromParameter(const Parameter &parameter);
mc::KineticMcSecondOmpi BuildKineticMcSecondOmpiFromParameter(const Parameter &parameter);
mc::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter);
mc::CanonicalMcSerial BuildCanonicalMcSerialFromParameter(const Parameter &parameter);
mc::CanonicalMcOmp BuildCanonicalMcOmpFromParameter(const Parameter &parameter);
ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter);
} // api

#endif //LMC_LMC_API_INCLUDE_HOME_H_

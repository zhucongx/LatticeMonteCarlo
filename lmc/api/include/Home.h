/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/6/23 3:05 PM                                                            *
 **************************************************************************************************/

#ifndef LMC_API_INCLUDE_HOME_H_
#define LMC_API_INCLUDE_HOME_H_
#include "Parameter.h"
#include "KineticMcFirstMpi.h"
#include "KineticMcFirstOmp.h"
#include "KineticMcChainOmpi.h"
#include "SimulatedAnnealing.h"
#include "CanonicalMcSerial.h"
#include "CanonicalMcOmp.h"
#include "Traverse.h"

namespace api {
void Print(const Parameter &parameter);
void Run(const Parameter &parameter);
mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter(const Parameter &parameter);
mc::KineticMcFirstOmp BuildKineticMcFirstOmpFromParameter(const Parameter &parameter);
mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter &parameter);
ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter);
mc::CanonicalMcSerial BuildCanonicalMcSerialFromParameter(const Parameter &parameter);
mc::CanonicalMcOmp BuildCanonicalMcOmpFromParameter(const Parameter &parameter);
ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter);
} // api

#endif //LMC_API_INCLUDE_HOME_H_

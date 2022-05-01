#ifndef LKMC_LKMC_API_INCLUDE_HOME_H_
#define LKMC_LKMC_API_INCLUDE_HOME_H_
#include "Parameter.h"
#include "ChainKmcMpi.h"
namespace api {
kmc::ChainKmcMpi BuildKmcMpiFromParameter(const Parameter &parameter);

} // namespace api

#endif //LKMC_LKMC_API_INCLUDE_HOME_H_

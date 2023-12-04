/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 6/16/23 1:53 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 5:20 PM                                                           *
 **************************************************************************************************/


/*! \file  Constants.hpp
 *  \brief File for the constants namespace definition.
 */

#ifndef LMC_CONSTANT_INCLUDE_CONSTANTS_HPP_
#define LMC_CONSTANT_INCLUDE_CONSTANTS_HPP_

/*! \brief Namespace for all constants.
 */
namespace constants {

/*! \brief Threshold to compare the values between the float numbers
 */
constexpr double kEpsilon = 1e-8;

/*! \brief Boltzmann constant in eV/K.
 */
constexpr double kBoltzmann = 8.617333262145e-5;

/*! \brief Prefactor used in kMC in Hz.
 */
constexpr double kPrefactor = 1e13; // Hz

} // constants
#endif //LMC_CONSTANT_INCLUDE_CONSTANTS_HPP_

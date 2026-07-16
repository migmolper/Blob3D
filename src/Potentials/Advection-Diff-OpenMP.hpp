/**
 * @file Potentials/Advection-Diff-OpenMP.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2023-05-19
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef ADVEC_DIFF_POTENTIAL_HPP
#define ADVEC_DIFF_POTENTIAL_HPP

#include "Macros.hpp"
#include <Eigen/Dense>

/**
 * @brief Function devoted to create the DMD function context of a Mg-Hx system
 *
 * @return governing_equations
 */
governing_equations Advection_Diff_constructor();


#endif /* ADVEC_DIFF_POTENTIAL_HPP */
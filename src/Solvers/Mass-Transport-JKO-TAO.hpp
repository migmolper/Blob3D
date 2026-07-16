/**
 * @file Kinetics/Mass-Transport-JKO-TAO.hpp
 * @author M.Molinos, MP.Ariza ([migmolper](https://github.com/migmolper) and
 * [mpariza](https://github.com/mpariza))
 * @brief
 * @version 0.1
 * @date 2026-03-17
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef MASS_TRANSPORT_JKO_TAO_HPP
#define MASS_TRANSPORT_JKO_TAO_HPP

#include "Macros.hpp"
#include <Eigen/Dense>

PetscErrorCode Mass_Trasport_Advection_Diffusion(
    PetscReal dt, DMD* Simulation, GoverningEquations& system_equations);

#endif /* MASS_TRANSPORT_JKO_TAO_HPP */

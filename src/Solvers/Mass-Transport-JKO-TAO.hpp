/**
 * @file Kinetics/Mass-Transport-JKO-TAO.hpp
 */

#ifndef MASS_TRANSPORT_JKO_TAO_HPP
#define MASS_TRANSPORT_JKO_TAO_HPP

#include "Simulation.hpp"
#include <Eigen/Dense>

PetscErrorCode
Mass_Trasport_Advection_Diffusion(PetscReal dt, Simulation &simulation,
                                  GoverningEquations &system_equations,
                                  boundaryCondition &boundary_conditions);

#endif /* MASS_TRANSPORT_JKO_TAO_HPP */

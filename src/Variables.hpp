/**
 * @file variables.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include "Macros.hpp"
#include <Eigen/Dense>

char OutputFolder[MAXC];
char InputFolder[MAXC];

/***********************************/
/********* PARALLELIZATION *********/
/***********************************/
PetscMPIInt rank_MPI = 0;
PetscMPIInt size_MPI = 1;

//! Number of MPI partitions in x,y,z directions
PetscInt ndiv_mesh_X = PETSC_DECIDE;
PetscInt ndiv_mesh_Y = PETSC_DECIDE;
PetscInt ndiv_mesh_Z = PETSC_DECIDE;

PetscInt size_MPI_X = PETSC_DECIDE; 
PetscInt size_MPI_Y = PETSC_DECIDE; 
PetscInt size_MPI_Z = PETSC_DECIDE;

/*********************************/
/** PARAMETERS FOR PETSC SOLVER **/
/*********************************/

//!< absolute convergence tolerance
double petsc_abstol = 1.e-7;

//!< relative convergence tolerance
double petsc_rtol = 1.e-10;

//!< convergence tolerance in terms of the norm of the change in the
//!< solution between steps, || delta x || < stol*|| x ||
double petsc_stol = 1.e-4;

//!< maximum number of iterations
double petsc_maxit = 30;

//!< maximum number of function evaluations
double petsc_maxf = 2000;

//!< Number of stored previous solutions and residuals
char petsc_ngmres_m[] = "2";

//!< The minimum step length
char petsc_linesearch_minlambda[] = "0.01";

//!< The linesearch damping parameter
char petsc_linesearch_damping[] = "0.50";

//!< The number of iterations for iterative line searches
char petsc_linesearch_max_it[] = "10";

/*********************************/
/*** PARAMETERS FOR DIFFUSION ***/
/*********************************/

//!< Diffusivity. Parameter to be adjusted. see article below
double Df = 0.045;

//!< time step. Parameter to be adjusted.
double dt_diffusion = 0.0004;

//!< diffusivity distance. Parameter to be adjusted.
double dr_diffusion = 3.2;

//!< label to consider diffusion or not. 0=No diffusion and
//!< 1= yes diffusion; by default yes diffusion
int diffusion = 1;

//!< minimun diffusion. If the diffusion between site i and j is less
//!< than this amoun, we dont take into account in our calculations.
double min_dxij = 1e-8;

//<! maximum number of iterations of the diffusion for each step
int max_it_diff = 1;

//<! in percentage
double max_total_mass = 0.01;
double max_Xh = 0.9;

/*****************************************/
DiffusivePotential Potential_AD;

/******************************************************/
/***** Variable to store the mass of each element *****/
/******************************************************/
double element_mass[112];

#endif /* VARIABLES_HPP */

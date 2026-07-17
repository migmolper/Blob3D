/**
 * @file variables.hpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief Global simulation parameters (declared here, defined in Variables.cpp).
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

extern char OutputFolder[MAXC];
extern char InputFolder[MAXC];

/***********************************/
/********* PARALLELIZATION *********/
/***********************************/
extern PetscMPIInt rank_MPI;
extern PetscMPIInt size_MPI;

//! Number of MPI partitions in x,y,z directions
extern PetscInt ndiv_mesh_X;
extern PetscInt ndiv_mesh_Y;
extern PetscInt ndiv_mesh_Z;

extern PetscInt size_MPI_X;
extern PetscInt size_MPI_Y;
extern PetscInt size_MPI_Z;

/*********************************/
/** PARAMETERS FOR PETSC SOLVER **/
/*********************************/

//!< absolute convergence tolerance
extern double petsc_abstol;

//!< relative convergence tolerance
extern double petsc_rtol;

//!< convergence tolerance in terms of the norm of the change in the
//!< solution between steps, || delta x || < stol*|| x ||
extern double petsc_stol;

//!< maximum number of iterations
extern double petsc_maxit;

//!< maximum number of function evaluations
extern double petsc_maxf;

//!< Number of stored previous solutions and residuals
extern char petsc_ngmres_m[];

//!< The minimum step length
extern char petsc_linesearch_minlambda[];

//!< The linesearch damping parameter
extern char petsc_linesearch_damping[];

//!< The number of iterations for iterative line searches
extern char petsc_linesearch_max_it[];

/*********************************/
/*** PARAMETERS FOR DIFFUSION ***/
/*********************************/

//!< Diffusivity. Parameter to be adjusted. see article below
extern double Df;

//!< time step. Parameter to be adjusted.
extern double dt_diffusion;

//!< diffusivity distance. Parameter to be adjusted.
extern double dr_diffusion;

//!< label to consider diffusion or not. 0=No diffusion and
//!< 1= yes diffusion; by default yes diffusion
extern int diffusion;

//!< minimun diffusion. If the diffusion between site i and j is less
//!< than this amoun, we dont take into account in our calculations.
extern double min_dxij;

//<! maximum number of iterations of the diffusion for each step
extern int max_it_diff;

//<! in percentage
extern double max_total_mass;
extern double max_Xh;



#endif /* VARIABLES_HPP */

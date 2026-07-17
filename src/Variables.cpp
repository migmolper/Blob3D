/**
 * @file Variables.cpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief Definitions of global simulation parameters.
 * @version 0.1
 * @date 2026-07-16
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "Variables.hpp"

char OutputFolder[MAXC];
char InputFolder[MAXC];

PetscMPIInt rank_MPI = 0;
PetscMPIInt size_MPI = 1;

PetscInt ndiv_mesh_X = PETSC_DECIDE;
PetscInt ndiv_mesh_Y = PETSC_DECIDE;
PetscInt ndiv_mesh_Z = PETSC_DECIDE;

PetscInt size_MPI_X = PETSC_DECIDE;
PetscInt size_MPI_Y = PETSC_DECIDE;
PetscInt size_MPI_Z = PETSC_DECIDE;

double petsc_abstol = 1.e-7;
double petsc_rtol = 1.e-10;
double petsc_stol = 1.e-4;
double petsc_maxit = 30;
double petsc_maxf = 2000;

char petsc_ngmres_m[] = "2";
char petsc_linesearch_minlambda[] = "0.01";
char petsc_linesearch_damping[] = "0.50";
char petsc_linesearch_max_it[] = "10";

double Df = 0.045;
double dt_diffusion = 0.0004;
double dr_diffusion = 3.2;
int diffusion = 1;
double min_dxij = 1e-8;
int max_it_diff = 1;
double max_total_mass = 0.01;
double max_Xh = 0.9;


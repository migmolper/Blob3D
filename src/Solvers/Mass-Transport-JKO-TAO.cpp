
/**
 * @file Kinetics/Mass-Transport-JKO-TAO.cpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief
 * @version 0.1
 * @date 2024-02-08
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <cstdlib>
#include <petscsystypes.h>
#if __APPLE__
#include <malloc/_malloc.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "Blobs/Blobs.hpp"
#include "Blobs/Ghosts.hpp"
#include "Blobs/Neighbors.hpp"
#include "Macros.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "petscdm.h"
#include "petscdmda.h"
#include "petscdmlabel.h"
#include <ctime>
#include <fstream>
#include <iomanip> // to print more decimals
#include <iostream>
#include <math.h>
#include <petscksp.h>
#include <petsctao.h>
#include <stdio.h>
#include <stdlib.h>

extern double petsc_abstol;
extern double petsc_rtol;
extern double petsc_stol;
extern double petsc_maxit;
extern double petsc_maxf;

struct JKO_ctx {

  /** @param blob_topology: List of neighs */
  const ParticleTopology *blob_topology;

  /*! @param box_idx_ptr: */
  const PetscInt *box_idx_ptr;

  /*! @param X_k: Reference value standard desviation of the
   * position */
  Vec X_k;

  /*! @param rho_k1: Density */
  Vec rho_k1;

  /*! @param beta_k1: Thermal lagrange multiplier*/
  Vec beta_k1;

  /*! @param beta_k: Thermal lagrange multiplier*/
  Vec beta_k;

  /*! @param mass: Mass vector */
  Vec mass;

  /*! @param system_equations Definition of the equation (non-owning) */
  GoverningEquations *system_equations;

  /*! @param boundary_conditions: Boundary conditions */
  boundaryCondition *boundary_conditions;

  /*! @param Delta_t: Time-step */
  PetscScalar Delta_t;

  /*! @param background_mesh:  */
  DM background_mesh;
};

static PetscErrorCode Advection(PetscReal dt, Simulation &simulation,
                                GoverningEquations &system_equations);

static PetscErrorCode JKO_Diffusion(PetscReal dt, Simulation &simulation,
                                    GoverningEquations &system_equations,
                                    boundaryCondition &boundary_conditions);

static PetscErrorCode compute_F0_and_RHS(Tao tao, Vec X_k1, PetscReal *F0,
                                         Vec D_JKO_Dx, void *ctx);

static PetscErrorCode compute_RHS(Tao tao, Vec X_k1, Vec D_JKO_Dx, void *ctx);

/************************************************************************/

PetscErrorCode
Mass_Trasport_Advection_Diffusion(PetscReal dt, Simulation &simulation,
                                  GoverningEquations &system_equations,
                                  boundaryCondition &boundary_conditions) {

  PetscFunctionBeginUser;

  PetscScalar Delta_r = 6.0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Advection step: update the particle positions according to U
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(Advection(dt, simulation, system_equations));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Source step: adjust the number of particles to account for sources
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(simulation.regenerate_topology(Delta_r));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Unconstrained diffusion step: update the particle position
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(
      JKO_Diffusion(dt, simulation, system_equations, boundary_conditions));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

static PetscErrorCode Advection(PetscReal dt, Simulation &simulation,
                                GoverningEquations &system_equations) {

  PetscFunctionBeginUser;
  (void)system_equations;
  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get system information
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! Get local number of sites in the simulation (without ghost)
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get local number of sites in the simulation (with ghost)
  PetscInt n_sites_local_ghosted;
  PetscCall(DMSwarmGetLocalSize(simulation.dm(), &n_sites_local_ghosted));

  //! Get number of ghost particles
  PetscInt n_sites_ghost = n_sites_local_ghosted - n_sites_local;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get list of mechanical neighbors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  IS *mechanical_neighs_idx = simulation.mechanical_neighs_idx();
  ParticleTopology *blob_topology = (ParticleTopology *)malloc(
      n_sites_local_ghosted * sizeof(ParticleTopology));

  for (PetscInt site_u = 0; site_u < n_sites_local_ghosted; site_u++) {
    PetscCall(DMSwarmGetParticleNeighbors(&blob_topology[site_u],
                                          mechanical_neighs_idx[site_u]));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Periodic box index
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt *box_idx_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "box-idx", NULL, NULL,
                            (void **)&box_idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index of the particles
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt *idx_q_ptr;
  PetscCall(
      DMSwarmGetField(simulation.dm(), "idx", NULL, NULL, (void **)&idx_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index for the ghost particles
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof_local = dim * n_sites_local;
  PetscInt n_dof_ghost = dim * n_sites_ghost;
  PetscInt *idx_dof_ghost = (PetscInt *)malloc(n_dof_ghost * sizeof(PetscInt));
  for (int i = 0; i < n_sites_ghost; i++) {
    for (int j = 0; j < dim; j++) {
      idx_dof_ghost[i * dim + j] = idx_q_ptr[n_sites_local + i] * dim + j;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get finite element mesh
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DM background_mesh;
  PetscCall(DMSwarmGetCellDM(simulation.dm(), &background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get mean position vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *mean_q_ptr; //! Mean position pointer
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void **)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get momentum vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *mean_p_ptr; //! Mean momentum pointer
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-p", NULL, NULL,
                            (void **)&mean_p_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get mass
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *mass_ptr; //! Mass pointer
  PetscCall(
      DMSwarmGetField(simulation.dm(), "mass", NULL, NULL, (void **)&mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Compute advection
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = 0; site_i < n_sites_local; site_i++) {
    for (int alpha = 0; alpha < dim; alpha++) {
      mean_p_ptr[site_i * dim + alpha] +=
          (mean_q_ptr[site_i * dim + alpha] / mass_ptr[site_i]) * dt;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Enforce periodic bcc and restore mean-q data
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec X; //! Mean position PETSc vector
  PetscCall(VecCreateGhostWithArray(PETSC_COMM_WORLD, n_dof_local,
                                    PETSC_DETERMINE, n_dof_ghost, idx_dof_ghost,
                                    mean_q_ptr, &X));
  PetscCall(VecSetUp(X));

  //! Update ghost values in the mean_q array
  PetscCall(VecGhostUpdateBegin(X, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(X, INSERT_VALUES, SCATTER_FORWARD));

  //! Enforce periodic bcc on the ghost particles
  PetscCall(VecEnforceGhostBlobsPeriodic(X, box_idx_ptr, background_mesh));

  //! Copy the updated ghost values back to the mean_q array
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL, NULL,
                                (void **)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore mean-p data
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-p", NULL, NULL,
                                (void **)&mean_p_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Restore mass data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mass", NULL, NULL,
                                (void **)&mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore Periodic box index
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "box-idx", NULL, NULL,
                                (void **)&box_idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore idx data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx", NULL, NULL,
                                (void **)&idx_q_ptr));
  free(idx_dof_ghost);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Free work space.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecDestroy(&X));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

static PetscErrorCode JKO_Diffusion(PetscReal dt, Simulation &simulation,
                                    GoverningEquations &system_equations,
                                    boundaryCondition &boundary_conditions) {

  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get system topology
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! Get local number of sites in the simulation (without ghost)
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get local number of sites in the simulation (with ghost)
  PetscInt n_sites_local_ghosted;
  PetscCall(DMSwarmGetLocalSize(simulation.dm(), &n_sites_local_ghosted));

  //! Get number of ghost particles
  PetscInt n_sites_ghost = n_sites_local_ghosted - n_sites_local;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get list of mechanical neighbors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  IS *mechanical_neighs_idx = simulation.mechanical_neighs_idx();
  ParticleTopology *blob_topology = (ParticleTopology *)malloc(
      n_sites_local_ghosted * sizeof(ParticleTopology));

  for (PetscInt site_u = 0; site_u < n_sites_local_ghosted; site_u++) {
    PetscCall(DMSwarmGetParticleNeighbors(&blob_topology[site_u],
                                          mechanical_neighs_idx[site_u]));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Periodic box index
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt *box_idx_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "box-idx", NULL, NULL,
                            (void **)&box_idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index of the particles
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt *idx_q_ptr;
  PetscCall(
      DMSwarmGetField(simulation.dm(), "idx", NULL, NULL, (void **)&idx_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index for the ghost particles
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof_local = dim * n_sites_local;
  PetscInt n_dof_ghost = dim * n_sites_ghost;
  PetscInt *idx_dof_ghost = (PetscInt *)malloc(n_dof_ghost * sizeof(PetscInt));
  for (int i = 0; i < n_sites_ghost; i++) {
    for (int j = 0; j < dim; j++) {
      idx_dof_ghost[i * dim + j] = idx_q_ptr[n_sites_local + i] * dim + j;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get index for the ghost particles
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ISLocalToGlobalMapping mapping_dofs;
  PetscCall(ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, dim, n_dof_local,
                                         idx_q_ptr, PETSC_USE_POINTER,
                                         &mapping_dofs));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get finite element mesh
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DM background_mesh;
  PetscCall(DMSwarmGetCellDM(simulation.dm(), &background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get mean position vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *mean_q_ptr; //! Mean position pointer
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void **)&mean_q_ptr));
  Vec X_k1; //! Mean position PETSc vector
  PetscCall(VecCreateGhostWithArray(PETSC_COMM_WORLD, n_dof_local,
                                    PETSC_DETERMINE, n_dof_ghost, idx_dof_ghost,
                                    mean_q_ptr, &X_k1));
  PetscCall(VecSetUp(X_k1));

#ifdef DEBUG_MODE
  PetscCall(VecView(X_k1, PETSC_VIEWER_STDOUT_WORLD));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get reference solution vector X-mean-q = {mean-q}
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec X_k;
  PetscCall(VecDuplicate(X_k1, &X_k));
  PetscCall(VecCopy(X_k1, X_k));
  PetscCall(VecGhostUpdateBegin(X_k, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(X_k, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecEnforceGhostBlobsPeriodic(X_k, box_idx_ptr, background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Initialize Jacobian matrix
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if FD_JACOBIAN == 1
  Mat Jac;
  PetscCall(MatCreate(PETSC_COMM_WORLD, &Jac));
  PetscCall(MatSetType(Jac, MATAIJ));
  PetscCall(MatSetSizes(Jac, n_dof_local, n_dof_local, PETSC_DETERMINE,
                        PETSC_DETERMINE));
  PetscCall(MatSetOptionsPrefix(Jac, "minV_dq_Jac_"));
  PetscCall(MatSetFromOptions(Jac));
  PetscCall(MatSetLocalToGlobalMapping(Jac, mapping_dofs, mapping_dofs));
  PetscCall(MatSetUp(Jac));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the molar fraction vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *rho_k1_ptr; //! Molr fraction pointer
  PetscCall(DMSwarmGetField(simulation.dm(), "rho", NULL, NULL,
                            (void **)&rho_k1_ptr));
  Vec rho_k1; //! Molar fraction PETSc vector
  PetscCall(VecCreateGhostWithArray(
      PETSC_COMM_WORLD, n_sites_local, PETSC_DETERMINE, n_sites_ghost,
      &idx_q_ptr[n_sites_local], rho_k1_ptr, &rho_k1));
  PetscCall(VecSetUp(rho_k1));
#ifdef DEBUG_MODE
  PetscCall(VecView(rho_k1, PETSC_VIEWER_STDOUT_WORLD));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the thermal multiplier vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *beta_k1_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "beta", NULL, NULL,
                            (void **)&beta_k1_ptr));
  Vec beta_k1; //! Thermal multiplier PETSc vector
  PetscCall(VecCreateGhostWithArray(
      PETSC_COMM_WORLD, n_sites_local, PETSC_DETERMINE, n_sites_ghost,
      &idx_q_ptr[n_sites_local], beta_k1_ptr, &beta_k1));
  PetscCall(VecSetUp(beta_k1));
#ifdef DEBUG_MODE
  PetscCall(VecView(beta_k1, PETSC_VIEWER_STDOUT_WORLD));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get reference solution vector X-mean-q = {mean-q}
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec beta_k;
  PetscCall(VecDuplicate(beta_k1, &beta_k));
  PetscCall(VecCopy(beta_k1, beta_k));
  PetscCall(VecGhostUpdateBegin(beta_k, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(beta_k, INSERT_VALUES, SCATTER_FORWARD));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the mass vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *mass_ptr; //! Mass pointer
  PetscCall(
      DMSwarmGetField(simulation.dm(), "mass", NULL, NULL, (void **)&mass_ptr));
  Vec mass; //! Mass PETSc vector
  PetscCall(VecCreateGhostWithArray(
      PETSC_COMM_WORLD, n_sites_local, PETSC_DETERMINE, n_sites_ghost,
      &idx_q_ptr[n_sites_local], mass_ptr, &mass));
  PetscCall(VecSetUp(mass));
#ifdef DEBUG_MODE
  PetscCall(VecView(mass, PETSC_VIEWER_STDOUT_WORLD));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define PETSc variables for the mechanical equilibrium
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Tao tao_min_JKO; //! Nonlinear optimization solver context
  KSP ksp;         //! linear solver context
  PC pc;           //! preconditioner context
  PetscInt TAO_iterations;
  // absolute convergence tolerance
  PetscReal abstol = 1e-3;
  // relative convergence tolerance
  PetscReal rtol = 1.e-10;
  // convergence tolerance in terms of the norm of the change in
  // the solution between steps, || delta x || < stol*|| x ||
  PetscReal stol = 1.e-6;
  // maximum number of iterations
  PetscInt maxit = 20;
  // maximum number of function evaluations
  PetscInt maxf = petsc_maxf;
  // Set a reason for convergence/divergence of solver (TAO)
  TaoConvergedReason min_JKO_reason;
  const char *TAO_strreason;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system; save final solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec D_JKO_Dx;
  PetscCall(VecDuplicate(X_k1, &D_JKO_Dx));
  PetscCall(VecSetOptionsPrefix(D_JKO_Dx, "minJKO_dx_RHS_"));
  PetscCall(VecSet(D_JKO_Dx, 0.0));
  PetscCall(VecGhostUpdateBegin(D_JKO_Dx, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(D_JKO_Dx, INSERT_VALUES, SCATTER_FORWARD));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     User-defined context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  struct JKO_ctx ctx;
  ctx.blob_topology = blob_topology;
  ctx.box_idx_ptr = box_idx_ptr;
  ctx.X_k = X_k;
  ctx.rho_k1 = rho_k1;
  ctx.beta_k1 = beta_k1;
  ctx.beta_k = beta_k;
  ctx.mass = mass;
  ctx.system_equations = &system_equations;
  ctx.boundary_conditions = &boundary_conditions;
  ctx.Delta_t = dt;
  ctx.background_mesh = background_mesh;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set solver parameters
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TaoCreate(PETSC_COMM_WORLD, &tao_min_JKO));
  PetscCall(TaoSetType(tao_min_JKO, TAOCG));
  PetscCall(TaoSetObjectiveAndGradient(tao_min_JKO, D_JKO_Dx,
                                       compute_F0_and_RHS, &ctx));
  PetscCall(TaoSetMaximumIterations(tao_min_JKO, maxit));
  PetscCall(TaoSetTolerances(tao_min_JKO, abstol, rtol, stol));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set solution vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TaoSetSolution(tao_min_JKO, X_k1));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Customize the non-linear solver (TAO)
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! Customized flag for this solver
  PetscCall(TaoSetOptionsPrefix(tao_min_JKO, "minJKO_dx_"));
  PetscCall(TaoSetFromOptions(tao_min_JKO));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Funtion to evaluate the Hessian using FD
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if FD_JACOBIAN == 1
  /*
     Color the matrix, i.e. determine groups of columns that share no common
    rows. These columns in the Jacobian can all be computed simultaneously.
  */
  PetscCall(MatColoringCreate(J, &coloring));
  PetscCall(MatColoringSetType(coloring, MATCOLORINGJP));
  PetscCall(MatColoringSetFromOptions(coloring));
  PetscCall(MatColoringApply(coloring, &iscoloring));
  PetscCall(MatColoringDestroy(&coloring));
  /*
     Create the data structure that SNESComputeJacobianDefaultColor() uses
     to compute the actual Jacobians via finite differences.
  */
  PetscCall(MatFDColoringCreate(J, iscoloring, &fdcoloring));
  PetscCall(MatFDColoringSetFunction(
      fdcoloring, (PetscErrorCode (*)(void))compute_RHS, &ctx));
  PetscCall(MatFDColoringSetUp(J, iscoloring, fdcoloring));
  PetscCall(ISColoringDestroy(&iscoloring));
  /*
    Tell SNES to use the routine SNESComputeJacobianDefaultColor()
    to compute Jacobians.
  */
  PetscCall(
      SNESSetJacobian(snes, J, J, SNESComputeJacobianDefaultColor, fdcoloring));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solve non-linear equations
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TaoSolve(tao_min_JKO));
  PetscCall(TaoGetTotalIterationNumber(tao_min_JKO, &TAO_iterations));
  PetscCall(TaoGetConvergedReason(tao_min_JKO, &min_JKO_reason));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                        "Number of TAO iterations = %" PetscInt_FMT "\n",
                        TAO_iterations));
  if (min_JKO_reason > 0) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Converged due to %s \n",
                          TaoConvergedReasons[min_JKO_reason]));
  } else if (min_JKO_reason <= 0) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Did not converged due to %s \n",
                          TaoConvergedReasons[min_JKO_reason]));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Enforce periodic bcc and restore mean-q data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostUpdateBegin(X_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(X_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecEnforceGhostBlobsPeriodic(X_k1, box_idx_ptr, background_mesh));
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL, NULL,
                                (void **)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore molar fraction data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "rho", NULL, NULL,
                                (void **)&rho_k1_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore thermal Lagrange Multiplier (beta) data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "beta", NULL, NULL,
                                (void **)&beta_k1_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore mass
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mass", NULL, NULL,
                                (void **)&mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore Periodic box index
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "box-idx", NULL, NULL,
                                (void **)&box_idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore idx data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx", NULL, NULL,
                                (void **)&idx_q_ptr));
  free(idx_dof_ghost);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore particle topology
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (PetscInt site_u = 0; site_u < n_sites_local_ghosted; site_u++) {
    PetscCall(DMSwarmRestoreParticleNeighbors(&blob_topology[site_u],
                                              mechanical_neighs_idx[site_u]));
  }
  free(blob_topology);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Free work space.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecDestroy(&D_JKO_Dx));
  PetscCall(VecDestroy(&X_k1));
  PetscCall(VecDestroy(&X_k));
  PetscCall(VecDestroy(&rho_k1));
  PetscCall(VecDestroy(&beta_k1));
  PetscCall(VecDestroy(&beta_k));
#if FD_JACOBIAN == 1
  PetscCall(MatDestroy(&J));
  PetscCall(MatFDColoringDestroy(&fdcoloring));
#endif
  PetscCall(TaoDestroy(&tao_min_JKO));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

static PetscErrorCode compute_F0_and_RHS(Tao tao, Vec X_k1,
                                         PetscReal *JKO_system, Vec D_JKO_Dx,
                                         void *ctx) {

  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get user context
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! Get the list of neighbors of each site
  const ParticleTopology *blob_topology = ((JKO_ctx *)ctx)->blob_topology;

  //! Get the periodic box index
  const PetscInt *box_idx_ptr = ((JKO_ctx *)ctx)->box_idx_ptr;

  //!
  Vec X_k = ((JKO_ctx *)ctx)->X_k;

  //! Get the density
  Vec rho_k1 = ((JKO_ctx *)ctx)->rho_k1;

  //!
  Vec beta_k1 = ((JKO_ctx *)ctx)->beta_k1;

  //!
  Vec beta_k = ((JKO_ctx *)ctx)->beta_k;

  //!
  Vec mass = ((JKO_ctx *)ctx)->mass;

  //! Take structure with the dmd equations
  GoverningEquations &system_equations = *((JKO_ctx *)ctx)->system_equations;

  //! Take structure with the boundary conditions
  boundaryCondition &boundary_conditions =
      *((JKO_ctx *)ctx)->boundary_conditions;

  //! Time step
  PetscScalar Delta_t = ((JKO_ctx *)ctx)->Delta_t;

  //! Get FE mesh
  DM background_mesh = ((JKO_ctx *)ctx)->background_mesh;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Update solution vector X_k1 and enforce periodic bcc
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostUpdateBegin(X_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(X_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecEnforceGhostBlobsPeriodic(X_k1, box_idx_ptr, background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Compute density
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(system_equations.evaluate_meassure_JKO(rho_k1, X_k1, beta_k1, mass,
                                                   blob_topology));
  PetscCall(VecGhostUpdateBegin(rho_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(rho_k1, INSERT_VALUES, SCATTER_FORWARD));

#ifdef DEBUG_MODE
  PetscCall(VecView(rho_k1, PETSC_VIEWER_STDOUT_WORLD));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Compute internal energy: JKO and barrier potential terms
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  *JKO_system = 0.0;
  PetscCall(system_equations.evaluate_JKO(JKO_system, Delta_t, rho_k1, X_k1,
                                          X_k, beta_k1, beta_k, mass,
                                          blob_topology));
  PetscCall(boundary_conditions.add_barrier_potential(JKO_system, X_k1, mass));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set the RHS to {mean_q_ref - mean_q}. This allow us to enforce the
    position for non-active sites. Then, loop in the subdomain of i_star
    to update RHS vector Y = {DV_Dmeanq} and enforce boundary conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! Evaluate the RHS of the JKO system
  PetscCall(system_equations.evaluate_D_JKO_Dq(
      D_JKO_Dx, Delta_t, rho_k1, X_k1, X_k, beta_k1, mass, blob_topology));
  //! Add the barrier forces to the JKO system
  PetscCall(boundary_conditions.add_barrier_forces(D_JKO_Dx, X_k1, mass));

  //! Update the ghost values of the RHS vector
  PetscCall(VecGhostUpdateBegin(D_JKO_Dx, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(D_JKO_Dx, INSERT_VALUES, SCATTER_FORWARD));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Output vectors if need it
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef DEBUG_MODE
  PetscCall(VecView(D_JKO_Dx, PETSC_VIEWER_STDOUT_WORLD));
#endif

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode compute_RHS(Tao tao, Vec X_k1, Vec D_JKO_Dq, void *ctx) {

  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get user context
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! Get the list of neighbors of each site
  const ParticleTopology *blob_topology = ((JKO_ctx *)ctx)->blob_topology;

  //! Get the periodic box index
  const PetscInt *box_idx_ptr = ((JKO_ctx *)ctx)->box_idx_ptr;

  //!
  Vec X_k = ((JKO_ctx *)ctx)->X_k;

  //! Get the density
  Vec rho_k1 = ((JKO_ctx *)ctx)->rho_k1;

  //!
  Vec beta_k1 = ((JKO_ctx *)ctx)->beta_k1;

  //!
  Vec mass = ((JKO_ctx *)ctx)->mass;

  //! Take structure with the dmd equations
  GoverningEquations &system_equations = *((JKO_ctx *)ctx)->system_equations;

  //! Take structure with the boundary conditions
  boundaryCondition &boundary_conditions =
      *((JKO_ctx *)ctx)->boundary_conditions;

  //! Time step
  PetscScalar Delta_t = ((JKO_ctx *)ctx)->Delta_t;

  //! Get FE mesh
  DM background_mesh = ((JKO_ctx *)ctx)->background_mesh;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Update solution vector X_k1 and enforce periodic bcc
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostUpdateBegin(X_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(X_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecEnforceGhostBlobsPeriodic(X_k1, box_idx_ptr, background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Compute density
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(system_equations.evaluate_meassure_JKO(rho_k1, X_k1, beta_k1, mass,
                                                   blob_topology));

  PetscCall(VecGhostUpdateBegin(rho_k1, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(rho_k1, INSERT_VALUES, SCATTER_FORWARD));

#ifdef DEBUG_MODE
  PetscCall(VecView(rho_k1, PETSC_VIEWER_STDOUT_WORLD));
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set the RHS to {mean_q_ref - mean_q}. This allow us to enforce the
    position for non-active sites. Then, loop in the subdomain of i_star
    to update RHS vector Y = {DV_Dmeanq} and enforce boundary conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(system_equations.evaluate_D_JKO_Dq(
      D_JKO_Dq, Delta_t, rho_k1, X_k1, X_k, beta_k1, mass, blob_topology));

  //! Add the barrier forces to the JKO system
  PetscCall(boundary_conditions.add_barrier_forces(D_JKO_Dq, X_k1, mass));

  PetscCall(VecGhostUpdateBegin(D_JKO_Dq, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(D_JKO_Dq, INSERT_VALUES, SCATTER_FORWARD));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Output vectors if need it
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef DEBUG_MODE
  PetscCall(VecView(D_JKO_Dq, PETSC_VIEWER_STDOUT_WORLD));
#endif

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/


/**
 * @file Mesh/test/test-Boundary-Conditions.cpp
 * @author Miguel Molinos (migmolper)
 * @brief
 * @version 0.1
 * @date 2025-02-04
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "Blobs/Blobs.hpp"
#include "Blobs/IO/dump-inputs.hpp"
#include "Blobs/IO/vtk-outputs.hpp"
#include "Blobs/Neighbors.hpp"
#include "Macros.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "Potentials/Advection-Diff-OpenMP.hpp"
#include "Solvers/Mass-Transport-JKO-TAO.hpp"
#include "Variables.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream> //std::cout//std::cin
#include <limits>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdlib.h>

extern int size_MPI;
extern int rank_MPI;

extern PetscInt ndiv_mesh_X;
extern PetscInt ndiv_mesh_Y;
extern PetscInt ndiv_mesh_Z;

extern PetscInt size_MPI_X;
extern PetscInt size_MPI_Y;
extern PetscInt size_MPI_Z;

//! Set help message
const char help[] = "BLOB3D blob particle method library (static)";

/********************************************************************************/

int main(int argc, char **argv) {

  unsigned int dim = NumberDimensions;
  const double Tolerance = 1.E-8;
  PetscScalar Delta_r = 6.0;
  PetscScalar kappa = 0.1;   //!
  PetscScalar Delta_t = 1;   //!
  PetscScalar m_p = 0.01;    //!
  PetscScalar rho_ref = 0.2; //!
  PetscErrorCode ierr;
  const char OutputFolder[MAXC] = "outputs";
  char SimulationFile[MAXC] = "inputs/Sphere-R-20.dump";
  char OutputFile[MAXC];

  ndiv_mesh_X = 60;
  ndiv_mesh_Y = 60;
  ndiv_mesh_Z = 60;

  size_MPI_X = 2;
  size_MPI_Y = 2;
  size_MPI_Z = 2;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  HPC libs
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//! Initialize MPI
#if USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
#endif

  //! Initialize PETSc/SLEPc
  PetscFunctionBeginUser;
#ifdef USE_SLEPC
  SlepcInitialize(&argc, &argv, (char *)0, help);
#else
  PetscInitialize(&argc, &argv, 0, help);
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Command line options
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_gatol", "1.e-6");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_gttol", "1.e-9");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_max_it", "100");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_type", "cg");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_cg_type", "prp");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_max_funcs", "100");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_monitor_globalization", "");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_view", "");
  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_converged_reason", "");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Read information from dump file
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  dump_file Simulation_dump_data = DMSwarmBlobsReadDump(SimulationFile);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Initialize atomistic simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Simulation simulation;
  PetscCall(simulation.initialize(
      Simulation_dump_data, BackgroundMeshType::DMDA_mesh, r_cutoff_default));
  PetscCall(
      DMSwarmSetMigrateType(simulation.dm(), DMSWARM_MIGRATE_DMCELLNSCATTER));
  PetscCall(DMSwarmMigrate(simulation.dm(), PETSC_TRUE));

  //! Set potential parameters
  Potential_AD.C = 0.0;
  Potential_AD.kappa = kappa;
  Potential_AD.rho_ref = rho_ref;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Free dump data
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMSwarmBlobsFreeDump(&Simulation_dump_data);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Output data
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMView(simulation.dm(), PETSC_VIEWER_STDOUT_WORLD));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set beta, mass and mean momentum value
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar beta_value = 2.0 / DSQR(Delta_r);
  Vec beta;
  PetscCall(DMSwarmCreateGlobalVectorFromField(simulation.dm(), "beta", &beta));
  PetscCall(VecSet(beta, beta_value));
  PetscCall(
      DMSwarmDestroyGlobalVectorFromField(simulation.dm(), "beta", &beta));

  Vec mass;
  PetscCall(DMSwarmCreateGlobalVectorFromField(simulation.dm(), "mass", &mass));
  PetscCall(VecSet(mass, m_p));
  PetscCall(
      DMSwarmDestroyGlobalVectorFromField(simulation.dm(), "mass", &mass));

  Vec P;
  PetscCall(DMSwarmCreateGlobalVectorFromField(simulation.dm(), "mean-p", &P));
  PetscCall(VecSet(P, 0.0));
  PetscCall(DMSwarmDestroyGlobalVectorFromField(simulation.dm(), "mean-p", &P));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create ghost atoms and topology
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(simulation.generate_topology(Delta_r));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Initialize  equations
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  AdvectionDiffusionEquations system_equations;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Outputs initial configuration
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  snprintf(OutputFile, MAXC, "%s/Mg-hcp-sphere-R-20-%i.vtp", OutputFolder, 0);
  PetscCall(DMSwarmBlobsViewVtk(simulation, OutputFile));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Update mean position
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (unsigned int i = 1; i < 50; i++) {

    //! Run advection-diffusion blob solver
    PetscCall(Mass_Trasport_Advection_Diffusion(Delta_t, simulation,
                                                system_equations));

    //! Update simulation topology
    PetscCall(simulation.regenerate_topology(Delta_r));

    //! Output
    snprintf(OutputFile, MAXC, "%s/Mg-hcp-sphere-R-20-%i.vtp", OutputFolder, i);
    PetscCall(DMSwarmBlobsViewVtk(simulation, OutputFile));
  }

  //! Destroy list of neighbors and other important information
  PetscCall(simulation.destroy_topology());

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Finalice HPC libs
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//! Finalize PETSc/SLEPs
#ifdef USE_SLEPC
  SlepcFinalize();
#else
  PetscFinalize();
#endif

//! Finalize MPI
#ifdef USE_MPI
  MPI_Finalize();
#endif
}

/********************************************************************************/
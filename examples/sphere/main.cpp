
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
#include "Boundaries/sphere.hpp"
#include "Macros.hpp"
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

  PetscErrorCode ierr;
  const char OutputFolder[MAXC] = "outputs";
  char SimulationFile[MAXC] = "inputs/Sphere-R-20.dump";
  char OutputFile[MAXC];

  ndiv_mesh_X = 60;
  ndiv_mesh_Y = 60;
  ndiv_mesh_Z = 60;

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

  size_MPI_X = atoi(argv[1]);
  size_MPI_Y = atoi(argv[2]);
  size_MPI_Z = atoi(argv[3]);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Number of MPI divisions: %i %i %i\n",
                        size_MPI_X, size_MPI_Y, size_MPI_Z));

  PetscInt NumberSteps = atoi(argv[4]);
  PetscCall(
      PetscPrintf(PETSC_COMM_WORLD, "Number of steps: %i\n", NumberSteps));

  PetscScalar kappa = atof(argv[5]); //!
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "kappa: %f\n", kappa));

  //! Initialize simulation parameters
  PetscInt NumberOfBlobs = 1477;
  PetscScalar TotalMass = 100;
  PetscScalar radius_domain = 50.0;
  PetscScalar radius_t0 = 20.0;
  PetscScalar volume_domain =
      (4.0 / 3.0) * M_PI * radius_domain * radius_domain * radius_domain;
  PetscScalar volume_init =
      (4.0 / 3.0) * M_PI * radius_t0 * radius_t0 * radius_t0;
  PetscScalar beta_value = 0.5 / DSQR(radius_domain);
  PetscScalar Delta_r = 2 * radius_domain * pow(NumberOfBlobs, -1.0 / 3.0); //!
  PetscScalar Delta_t = Delta_r * Delta_r / kappa;                          //!
  PetscScalar penalty = 1.0 / Delta_t;
  PetscScalar penalty_buffer = 2.0;
  PetscScalar m_p = TotalMass / NumberOfBlobs;     //!
  PetscScalar rho_ref = TotalMass / volume_domain; //!
  PetscScalar rho_init = TotalMass / volume_init;  //!

  if (rank_MPI == 0) {
    std::cout << "Delta_r: " << Delta_r << std::endl;
    std::cout << "Delta_t: " << Delta_t << std::endl;
    std::cout << "penalty: " << penalty << std::endl;
    std::cout << "m_p: " << m_p << std::endl;
    std::cout << "rho_ref: " << rho_ref << std::endl;
    std::cout << "rho_init: " << rho_init << std::endl;
    std::cout << "beta_value: " << beta_value << std::endl;
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Read information from dump file
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  dump_file Simulation_dump_data = DMSwarmBlobsReadDump(SimulationFile);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Initialize blob simulation
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Simulation simulation;
  PetscCall(simulation.initialize(Simulation_dump_data,
                                  BackgroundMeshType::DMDA_mesh, Delta_r));
  PetscCall(
      DMSwarmSetMigrateType(simulation.dm(), DMSWARM_MIGRATE_DMCELLNSCATTER));
  PetscCall(DMSwarmMigrate(simulation.dm(), PETSC_TRUE));

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
  Vec beta;
  PetscCall(DMSwarmCreateGlobalVectorFromField(simulation.dm(), "beta", &beta));
  PetscCall(VecSet(beta, beta_value));
  PetscCall(
      DMSwarmDestroyGlobalVectorFromField(simulation.dm(), "beta", &beta));

  Vec rho;
  PetscCall(DMSwarmCreateGlobalVectorFromField(simulation.dm(), "rho", &rho));
  PetscCall(VecSet(rho, rho_init));
  PetscCall(DMSwarmDestroyGlobalVectorFromField(simulation.dm(), "rho", &rho));

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
   Create ghost blobs and topology
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(simulation.generate_topology(Delta_r));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Initialize  equations
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  AdvectionDiffusionEquations system_equations(kappa, rho_ref);
  BC_Sphere boundary_conditions(radius_domain, Eigen::Vector3d(0.0, 0.0, 0.0),
                                penalty, penalty_buffer);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Outputs initial configuration
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  snprintf(OutputFile, MAXC, "%s/Mg-hcp-sphere-R-20-%i.vtp", OutputFolder, 0);
  PetscCall(DMSwarmBlobsViewVtk(simulation, OutputFile));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Update mean position
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (unsigned int i = 1; i < NumberSteps; i++) {

    //! Run advection-diffusion blob solver
    PetscCall(Mass_Trasport_Advection_Diffusion(
        Delta_t, simulation, system_equations, boundary_conditions));

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
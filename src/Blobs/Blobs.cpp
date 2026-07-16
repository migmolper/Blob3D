/**
 * @file Atoms/Atom.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <petscis.h>
#include <petscsystypes.h>
#if __APPLE__
#include <malloc/_malloc.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "Atoms/Atom.hpp"
#include "Macros.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "Mesh/Create-Mesh.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern PetscMPIInt size_MPI;
extern PetscMPIInt rank_MPI;

extern PetscInt ndiv_mesh_X;
extern PetscInt ndiv_mesh_Y;
extern PetscInt ndiv_mesh_Z;

using namespace std;

/*******************************************************/

static PetscErrorCode InitializeBlobs_DMDA(DM* atomistic_data,
                                           dump_file Simulation_file) {

  PetscFunctionBeginUser;

  //! Get Integers
  unsigned int dim = NumberDimensions;
  PetscInt n_atoms = Simulation_file.n_atoms;
  PetscInt n_atoms_local = 0;
  MPI_Comm comm;
  PetscMPIInt size;
  PetscMPIInt rank;

  PetscCall(PetscObjectGetComm((PetscObject)*atomistic_data, &comm));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCallMPI(MPI_Comm_rank(comm, &rank));

  //! @brief Get the bounding box of the background mesh
  DM FE_Mesh;
  PetscReal ll_vertex[3], ur_vertex[3];
  PetscCall(DMSwarmGetCellDM(*atomistic_data, &FE_Mesh));
  PetscCall(DMGetLocalBoundingBox(FE_Mesh, ll_vertex, ur_vertex));

  //! @brief mean_q: Mean value of each atomic position
  double* mean_q_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, DMSwarmPICField_coor, NULL, NULL,
                            (void**)&mean_q_ptr));
  Eigen::Map<MatrixType> mean_q(mean_q_ptr, n_atoms, dim);

  //! @brief idx_ptr: Global index of each atomic position
  int* idx_ptr;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  //! @brief Loop in the elemnts of the background mesh and find atoms in the
  //! mesh
  for (PetscInt site_i = 0; site_i < n_atoms; site_i++) {

    Eigen::Map<const Eigen::Vector3d> mean_q_i(
        &(Simulation_file.mean_q[site_i * dim]), dim);

    if ((mean_q_i(0) >= ll_vertex[0]) &&  //!
        (mean_q_i(0) < ur_vertex[0]) &&   //!
        (mean_q_i(1) >= ll_vertex[1]) &&  //!
        (mean_q_i(1) < ur_vertex[1]) &&   //!
        (mean_q_i(2) >= ll_vertex[2]) &&  //!
        (mean_q_i(2) < ur_vertex[2])) {
      idx_ptr[n_atoms_local] = site_i;
      mean_q.row(n_atoms_local) = mean_q_i;
      n_atoms_local++;
    }
  }

  PetscCall(DMSwarmRestoreField(*atomistic_data, DMSwarmPICField_coor, NULL,
                                NULL, (void**)&mean_q_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "idx", NULL, NULL,
                                (void**)&idx_ptr));

  //! Check global size
  PetscInt n_atoms_all = 0;
  PetscCall(MPIU_Allreduce(&n_atoms_local, &n_atoms_all, 1, MPI_INT, MPI_SUM,
                           MPI_COMM_WORLD));

  if (n_atoms_all != n_atoms) {
    PetscCall(PetscError(
        PETSC_COMM_WORLD, __LINE__, "init_Blob_simulation", __FILE__,
        PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
        "The global size (%i) of the DMSwarm does not match the input (%i) !",
        n_atoms_all, n_atoms));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  //! @brief Set local sizes
  PetscCall(DMSwarmSetLocalSizes(*atomistic_data, n_atoms_local, BufferLenght));

  //! @brief mean_q_ref: Reference value of the mean position of each atomic
  //! site
  double* mean_q_ref;
  PetscCall(DMSwarmGetField(*atomistic_data, "mean-q-ref", NULL, NULL,
                            (void**)&mean_q_ref));

  //! @brief mean_p: Mean value of each atomic momentum
  double* mean_p;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "mean-p", NULL, NULL, (void**)&mean_p));

  //! @brief rho: Density
  double* rho;
  PetscCall(DMSwarmGetField(*atomistic_data, "rho", NULL, NULL, (void**)&rho));

  //! @brief mass: Mass of each particle
  double* mass;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "mass", NULL, NULL, (void**)&mass));

  //! @brief beta: Thermal Lagrange multiplier
  double* beta;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "beta", NULL, NULL, (void**)&beta));

  //! @brief idx_ptr: index of the site
  PetscCall(
      DMSwarmGetField(*atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  //! @brief ghost_ptr: index if the site is a ghost atom or not
  PetscInt* ghost_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, "ghost", NULL, NULL,
                            (void**)&ghost_ptr));

  //! @brief box_idx_ptr: periodic box index
  PetscInt* box_idx_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, "box-idx", NULL, NULL,
                            (void**)&box_idx_ptr));

  //! @brief site_mpi_rank: Integer which defines MPI rank location
  PetscInt* site_mpi_rank;
  PetscCall(DMSwarmGetField(*atomistic_data, "MPI-rank", NULL, NULL,
                            (void**)&site_mpi_rank));

  //!  @brief Pointers with the information of the boundary condition
  PetscInt* idx_bcc_mean_q_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, "idx-bcc-mean-q", NULL, NULL,
                            (void**)&idx_bcc_mean_q_ptr));

  PetscInt* idx_bcc_beta_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, "idx-bcc-beta", NULL, NULL,
                            (void**)&idx_bcc_beta_ptr));

  //! Set information from the .dump file
  for (PetscInt local_site_i = 0; local_site_i < n_atoms_local;
       local_site_i++) {

    PetscInt site_i = idx_ptr[local_site_i];

    site_mpi_rank[local_site_i] = (PetscInt)rank;

    ghost_ptr[local_site_i] = 0;

    box_idx_ptr[local_site_i] = 0;

    for (int alpha = 0; alpha < dim; alpha++) {
      mean_q_ref[local_site_i * dim + alpha] =
          Simulation_file.mean_q[site_i * dim + alpha];
      mean_p[local_site_i * dim + alpha] = 0.0;
    }

    rho[local_site_i] = Simulation_file.xi[site_i];
    if ((rho[local_site_i] < 0) || (rho[local_site_i] > 1.0)) {
      PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "InitializeBlobs_DMDA",
                           __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                           "The density at site %i should take value "
                           "between 0 and 1. Value: %f",
                           site_i, rho[local_site_i]));
      PetscFunctionReturn(PETSC_ERR_RETURN);
    }

    beta[local_site_i] = Simulation_file.beta[site_i];
    if (beta[local_site_i] < 0) {
      PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "InitializeBlobs_DMDA",
                           __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                           "The thermal multiplier at site %i should be larger "
                           "than zero. Value: %f",
                           site_i, beta[local_site_i]));
      PetscFunctionReturn(PETSC_ERR_RETURN);
    }

    mass[local_site_i] = 0.0;

    idx_bcc_mean_q_ptr[local_site_i] = 0;

    idx_bcc_beta_ptr[local_site_i] = Simulation_file.beta_bcc[site_i];
    if ((idx_bcc_beta_ptr[local_site_i] != 0) &&
        (idx_bcc_beta_ptr[local_site_i] != 1)) {
      PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "InitializeBlobs_DMDA",
                           __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                           "The bcc indicator for the thermal multiplier at "
                           "site %i should be 0 or 1",
                           site_i));
      PetscFunctionReturn(PETSC_ERR_RETURN);
    }
  }

  //! Restore particle fields
  PetscCall(DMSwarmRestoreField(*atomistic_data, "mean-q-ref", NULL, NULL,
                                (void**)&mean_q_ref));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "mean-p", NULL, NULL,
                                (void**)&mean_p));
  PetscCall(
      DMSwarmRestoreField(*atomistic_data, "rho", NULL, NULL, (void**)&rho));
  PetscCall(
      DMSwarmRestoreField(*atomistic_data, "mass", NULL, NULL, (void**)&mass));
  PetscCall(
      DMSwarmRestoreField(*atomistic_data, "beta", NULL, NULL, (void**)&beta));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "idx", NULL, NULL,
                                (void**)&idx_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "ghost", NULL, NULL,
                                (void**)&ghost_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "box-idx", NULL, NULL,
                                (void**)&box_idx_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "MPI-rank", NULL, NULL,
                                (void**)&site_mpi_rank));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "idx-bcc-mean-q", NULL, NULL,
                                (void**)&idx_bcc_mean_q_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "idx-bcc-beta", NULL, NULL,
                                (void**)&idx_bcc_beta_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/*******************************************************/

static PetscErrorCode InitializeBlobs_DMPlex(DM* atomistic_data,
                                             dump_file Simulation_file) {

  PetscFunctionBeginUser;

  //! Get Integers
  unsigned int dim = NumberDimensions;
  PetscInt n_atoms = Simulation_file.n_atoms;
  PetscInt n_atoms_local = 0;
  MPI_Comm comm;
  PetscMPIInt size;
  PetscMPIInt rank;

  PetscCall(PetscObjectGetComm((PetscObject)*atomistic_data, &comm));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCallMPI(MPI_Comm_rank(comm, &rank));

  //! @brief Get the bounding box of the background mesh
  DM FE_Mesh;
  PetscReal ll_vertex[3], ur_vertex[3];
  PetscCall(DMSwarmGetCellDM(*atomistic_data, &FE_Mesh));
  PetscCall(DMGetLocalBoundingBox(FE_Mesh, ll_vertex, ur_vertex));

  //! @brief Check local dimensions (to-do)

  //! @brief mean_q: Mean value of each atomic position
  double* mean_q_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, DMSwarmPICField_coor, NULL, NULL,
                            (void**)&mean_q_ptr));
  Eigen::Map<MatrixType> mean_q(mean_q_ptr, n_atoms, dim);

  //! @brief idx_ptr: Global index of each atomic position
  int* idx_ptr;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  //! @brief Loop in the elemnts of the background mesh and find atoms in the
  //! mesh
  for (PetscInt site_i = 0; site_i < n_atoms; site_i++) {

    Eigen::Map<const Eigen::Vector3d> mean_q_i(
        &(Simulation_file.mean_q[site_i * dim]), dim);

    if ((mean_q_i(0) >= ll_vertex[0]) &&  //!
        (mean_q_i(0) < ur_vertex[0]) &&   //!
        (mean_q_i(1) >= ll_vertex[1]) &&  //!
        (mean_q_i(1) < ur_vertex[1]) &&   //!
        (mean_q_i(2) >= ll_vertex[2]) &&  //!
        (mean_q_i(2) < ur_vertex[2])) {
      idx_ptr[n_atoms_local] = site_i;
      mean_q.row(n_atoms_local) = mean_q_i;
      n_atoms_local++;
    }
  }

  PetscCall(DMSwarmRestoreField(*atomistic_data, DMSwarmPICField_coor, NULL,
                                NULL, (void**)&mean_q_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "idx", NULL, NULL,
                                (void**)&idx_ptr));

  //! Check global size
  PetscInt n_atoms_all = 0;
  PetscCall(MPIU_Allreduce(&n_atoms_local, &n_atoms_all, 1, MPI_INT, MPI_SUM,
                           MPI_COMM_WORLD));

  if (n_atoms_all != n_atoms) {
    PetscCall(PetscError(
        PETSC_COMM_WORLD, __LINE__, "init_Blob_simulation", __FILE__,
        PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
        "The global size (%i) of the DMSwarm does not match the input (%i) !",
        n_atoms_all, n_atoms));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  //! @brief Set local sizes
  PetscCall(DMSwarmSetLocalSizes(*atomistic_data, n_atoms_local, BufferLenght));

  //! @brief mean_q_ref: Reference value of the mean position of each atomic
  //! site
  double* mean_q_ref;
  PetscCall(DMSwarmGetField(*atomistic_data, "mean-q-ref", NULL, NULL,
                            (void**)&mean_q_ref));

  //! @brief mean_p: Mean value of each atomic momentum
  double* mean_p;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "mean-p", NULL, NULL, (void**)&mean_p));

  //! @brief rho: Density
  double* rho;
  PetscCall(DMSwarmGetField(*atomistic_data, "rho", NULL, NULL, (void**)&rho));

  //! @brief mass: Mass of each particle
  double* mass;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "mass", NULL, NULL, (void**)&mass));

  //! @brief beta: Thermal Lagrange multiplier
  double* beta;
  PetscCall(
      DMSwarmGetField(*atomistic_data, "beta", NULL, NULL, (void**)&beta));

  //! @brief idx_ptr: index of the site
  PetscCall(
      DMSwarmGetField(*atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  //! @brief ghost_ptr: index if the site is a ghost atom or not
  PetscInt* ghost_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, "ghost", NULL, NULL,
                            (void**)&ghost_ptr));

  //! @brief box_idx_ptr: periodic box index
  PetscInt* box_idx_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, "box-idx", NULL, NULL,
                            (void**)&box_idx_ptr));

  //! @brief site_mpi_rank: Integer which defines MPI rank location
  PetscInt* site_mpi_rank;
  PetscCall(DMSwarmGetField(*atomistic_data, "MPI-rank", NULL, NULL,
                            (void**)&site_mpi_rank));

  //! @brief Pointers with the information of the boundary condition
  PetscInt* idx_bcc_beta_ptr;
  PetscCall(DMSwarmGetField(*atomistic_data, "idx-bcc-beta", NULL, NULL,
                            (void**)&idx_bcc_beta_ptr));

  //! Set information from the .dump file
  for (PetscInt local_site_i = 0; local_site_i < n_atoms_local;
       local_site_i++) {

    PetscInt site_i = idx_ptr[local_site_i];

    site_mpi_rank[local_site_i] = (PetscInt)rank;

    ghost_ptr[local_site_i] = 0;

    box_idx_ptr[local_site_i] = 0;

    for (int alpha = 0; alpha < dim; alpha++) {
      mean_q_ref[local_site_i * dim + alpha] =
          Simulation_file.mean_q[site_i * dim + alpha];
      mean_p[local_site_i * dim + alpha] = 0.0;
    }

    rho[local_site_i] = Simulation_file.xi[site_i];
    if ((rho[local_site_i] < 0) || (rho[local_site_i] > 1.0)) {
      PetscCall(PetscError(
          PETSC_COMM_SELF, __LINE__, "InitializeBlobs_DMPlex", __FILE__,
          PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
          "The density at site %i should take value between 0 and 1",
          local_site_i));
      PetscFunctionReturn(PETSC_ERR_RETURN);
    }

    mass[local_site_i] = 0.0;

    beta[local_site_i] = Simulation_file.beta[site_i];
    if (beta[local_site_i] < 0) {
      PetscCall(PetscError(
          PETSC_COMM_SELF, __LINE__, "InitializeBlobs_DMPlex", __FILE__,
          PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
          "The thermal multiplier at site %i should be larger than zero",
          local_site_i));
      PetscFunctionReturn(PETSC_ERR_RETURN);
    }

    idx_bcc_beta_ptr[local_site_i] = Simulation_file.beta_bcc[site_i];
    if ((idx_bcc_beta_ptr[local_site_i] != 0) ||
        (idx_bcc_beta_ptr[local_site_i] != 1)) {
      PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "InitializeBlobs_DMPlex",
                           __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                           "The bcc indicator for the thermal multiplier at "
                           "site %i should be 0 or 1",
                           local_site_i));
      PetscFunctionReturn(PETSC_ERR_RETURN);
    }
  }

  //! Restore particle fields
  PetscCall(DMSwarmRestoreField(*atomistic_data, "mean-q-ref", NULL, NULL,
                                (void**)&mean_q_ref));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "mean-p", NULL, NULL,
                                (void**)&mean_p));
  PetscCall(
      DMSwarmRestoreField(*atomistic_data, "rho", NULL, NULL, (void**)&rho));
  PetscCall(
      DMSwarmRestoreField(*atomistic_data, "mass", NULL, NULL, (void**)&mass));
  PetscCall(
      DMSwarmRestoreField(*atomistic_data, "beta", NULL, NULL, (void**)&beta));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "idx", NULL, NULL,
                                (void**)&idx_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "ghost", NULL, NULL,
                                (void**)&ghost_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "box-idx", NULL, NULL,
                                (void**)&box_idx_ptr));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "MPI-rank", NULL, NULL,
                                (void**)&site_mpi_rank));
  PetscCall(DMSwarmRestoreField(*atomistic_data, "idx-bcc-beta", NULL, NULL,
                                (void**)&idx_bcc_beta_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/*******************************************************/

/*
 Create a DMShell and attach a regularly spaced DMDA for point location
 Override methods for point location
*/
PetscErrorCode init_Blob_simulation(DMD* Simulation,               //!
                                    dump_file Simulation_file,     //!
                                    BackgroundMeshType mesh_type,  //!
                                    double r_cutoff_V) {

  unsigned int dim = NumberDimensions;

  DM atomistic_data;
  DM FE_Mesh;
  MPI_Comm comm;
  PetscMPIInt size;
  PetscMPIInt rank;
  PetscInt blocksize;

  //!
  PetscInt n_atoms_global = 0;
  PetscInt n_atoms_local = 0;
  PetscInt n_atoms = Simulation_file.n_atoms;

  //!
  IS beta_bcc;
  IS gamma_bcc;

  PetscFunctionBeginUser;

  /* Create the bg FE */
  PetscCall(CreateMesh(&FE_Mesh, Simulation_file, mesh_type, r_cutoff_V));

  /* Create the swarm */
  PetscCall(PetscObjectGetComm((PetscObject)FE_Mesh, &comm));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  PetscCall(DMCreate(comm, &atomistic_data));
  PetscCall(DMSetType(atomistic_data, DMSWARM));
  PetscCall(DMSetDimension(atomistic_data, dim));
  PetscCall(PetscObjectSetName((PetscObject)atomistic_data, "Atoms"));
  PetscCall(DMSwarmSetType(atomistic_data, DMSWARM_PIC));
  PetscCall(DMSwarmSetCellDM(atomistic_data, FE_Mesh));

  //! Molar fraction: Xi := <n>_0
  PetscCall(
      DMSwarmRegisterPetscDatatypeField(atomistic_data, "rho", 1, PETSC_REAL));

  //! Thermal multiplier: Beta := 1/(kb * T)
  PetscCall(
      DMSwarmRegisterPetscDatatypeField(atomistic_data, "beta", 1, PETSC_REAL));

  //! Mass of the particle
  PetscCall(
      DMSwarmRegisterPetscDatatypeField(atomistic_data, "mass", 1, PETSC_REAL));

  //! Refernece value of the mean value of the position: mean-q-ref
  PetscCall(DMSwarmRegisterPetscDatatypeField(atomistic_data, "mean-q-ref", 3,
                                              PETSC_REAL));

  //! Mean momentum of the particle
  PetscCall(DMSwarmRegisterPetscDatatypeField(atomistic_data, "mean-p", 3,
                                              PETSC_REAL));

  //! Site ghost
  PetscCall(
      DMSwarmRegisterPetscDatatypeField(atomistic_data, "ghost", 1, PETSC_INT));

  //! Periodic box index
  PetscCall(DMSwarmRegisterPetscDatatypeField(atomistic_data, "box-idx", 1,
                                              PETSC_INT));

  //! Site idx
  PetscCall(
      DMSwarmRegisterPetscDatatypeField(atomistic_data, "idx", 1, PETSC_INT));

  //! MPI rank
  PetscCall(DMSwarmRegisterPetscDatatypeField(atomistic_data, "MPI-rank", 1,
                                              PETSC_INT));

  //! Boundary condition idex for the mean position
  PetscCall(DMSwarmRegisterPetscDatatypeField(atomistic_data, "idx-bcc-mean-q",
                                              1, PETSC_INT));

  //! Boundary condition idex for the shape of the distribution
  PetscCall(DMSwarmRegisterPetscDatatypeField(atomistic_data, "idx-bcc-beta", 1,
                                              PETSC_INT));

  //!
  PetscCall(DMSwarmFinalizeFieldRegister(atomistic_data));
  PetscCall(DMSwarmSetLocalSizes(atomistic_data, n_atoms, BufferLenght));

  //! Initialize the particles using dump information
  PetscCall(InitializeBlobs_DMDA(&atomistic_data, Simulation_file));

  //! Set the number of particles
  PetscCall(DMSwarmGetSize(atomistic_data, &n_atoms_global));
  PetscCall(DMSwarmGetLocalSize(atomistic_data, &n_atoms_local));

  //! Create a context to relate petsc numbering and dump numbering
  AO dump2petsc_mapping;
  IS dump_ordering;

  PetscInt* idx_ptr;
  PetscCall(
      DMSwarmGetField(atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, n_atoms_local, idx_ptr,
                            PETSC_COPY_VALUES, &dump_ordering));
  PetscCall(AOCreateBasicIS(dump_ordering, NULL, &dump2petsc_mapping));
  PetscCall(ISDestroy(&dump_ordering));
  PetscCall(AOApplicationToPetsc(dump2petsc_mapping, n_atoms_local, idx_ptr));
  PetscCall(
      DMSwarmRestoreField(atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  //! Fill the context structure
  Simulation->n_sites_global = n_atoms_global;
  Simulation->n_sites_local = n_atoms_local;
  Simulation->Temperature_env = 0.0;
  Simulation->ChemicalPotential_env = 0.0;
  Simulation->Pressure_env = 0.0;
  Simulation->F = Eigen::Matrix3d::Identity();
  Simulation->atomistic_data = atomistic_data;
  Simulation->dump2petsc_mapping = dump2petsc_mapping;

  PetscFunctionReturn(PETSC_SUCCESS);
}

/*******************************************************/

PetscErrorCode destroy_Blob_simulation(DMD* Simulation) {

  PetscFunctionBeginUser;

  //! @brief Destroy the DMSwarm
  DM FE_Mesh;
  PetscCall(DMSwarmGetCellDM(Simulation->atomistic_data, &FE_Mesh));
  PetscCall(DMDestroy(&FE_Mesh));

  //! @brief Atomistic variables
  PetscCall(DMDestroy(&Simulation->atomistic_data));

  //! @brief Topological variables
  PetscCall(AODestroy(&Simulation->dump2petsc_mapping));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/*******************************************************/

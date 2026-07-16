/**
 * @file boundary_conditions.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Mesh/Boundary-Conditions.hpp"
#include "Blobs/Blobs.hpp"
#include "Blobs/Neighbors.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Macros.hpp"
#include "Mesh/Coordinates.hpp"
#include "Mesh/InOut-Mesh.hpp"
#include "petscsys.h"
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <petscsystypes.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/************************************************************************/

/**
 * @brief Fold a physical coordinate into the fundamental cell [lo, hi).
 *
 * Used before DMDA/DMSwarm partition tests so periodic image atoms match the
 * half-open brick owned by PETSc ranks.
 */
inline PetscReal FoldPeriodicCoordinate(PetscReal x, PetscReal lo,
                                        PetscReal hi) {
  const PetscReal L = hi - lo;
  if (L <= 0.0) return x;
  x -= lo;
  x -= L * floor(x / L);
  if (x >= L) x = 0.0;
  if (x < 0.0) x = 0.0;
  return lo + x;
}

/************************************************************************/

PetscErrorCode get_mesh_boundary_condition(DMBoundaryType* bcc, DM* da) {
  PetscFunctionBeginUser;

  DMBoundaryType bx, by, bz;
  PetscCall(DMDAGetInfo(*da, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                        NULL, &bx, &by, &bz, NULL));

  bcc[0] = bx;
  bcc[1] = by;
  bcc[2] = bz;

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmEnforceAtomsPeriodic(Simulation& simulation,
                                           double buffer_width) {

  PetscFunctionBeginUser;
  (void)buffer_width;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  unsigned int dim = NumberDimensions;

  //! Number of local physical sites
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get atomistic coordinates
  PetscScalar* mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the global (x,y,z) indices of the lower left corner and size of the
    local region, excluding ghost points.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscReal gmin[3], gmax[3];
  DM background_mesh;
  PetscCall(DMSwarmGetCellDM(simulation.dm(), &background_mesh));
  PetscCall(DMGetPhysicalBoundingBox(background_mesh, gmin, gmax, NULL));

  //! Define computational box
  const PetscScalar box_coords[6] = {gmin[0], gmin[1], gmin[2],
                                     gmax[0], gmax[1], gmax[2]};

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get mesh boundary conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMBoundaryType bcc[3] = {DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DM_BOUNDARY_NONE};
  PetscCall(get_mesh_boundary_condition(bcc, &background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check if the site is inside the domain and enforce periodic boundary
  conditions
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! @brief found: Flag to check if the atom is found in the mesh
  int* found_local_ptr = (int*)calloc(n_sites_local, sizeof(int));

#pragma omp parallel for schedule(runtime)
  for (PetscInt site_i = 0; site_i < n_sites_local; site_i++) {

    //! Update site position with the crystal direction
    Eigen::Vector3d mean_q_i;
    mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
        mean_q_ptr[site_i * dim + 2];

    //! Check if the updated site is inside the domain
    if (In_Out_Mesh_Closed(mean_q_i, box_coords) == true) {
      found_local_ptr[site_i] = 1;
    } else {
      //! Apply periodic boundary conditions
      for (unsigned int alpha = 0; alpha < dim; alpha++) {
        if (bcc[alpha] == DM_BOUNDARY_PERIODIC) {
          mean_q_ptr[site_i * dim + alpha] =
              FoldPeriodicCoordinate(mean_q_i(alpha), gmin[alpha], gmax[alpha]);
        } else {
          mean_q_ptr[site_i * dim + alpha] = mean_q_i(alpha);
        }
      }
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Check if the site is inside the domain after the enforcement of periodic
  boundary conditions if not, error is raised
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (PetscInt site_i = 0; site_i < n_sites_local; site_i++) {
    if (found_local_ptr[site_i] == 0) {

      //! Update site position with the crystal direction
      Eigen::Vector3d mean_q_i;
      mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
          mean_q_ptr[site_i * dim + 2];

      //! Check if the updated site is inside the domain
      if (In_Out_Mesh_Closed(mean_q_i, box_coords) == true) {
        found_local_ptr[site_i] = 1;
      } else {
        PetscCall(PetscError(
            PETSC_COMM_SELF, __LINE__, "DMSwarmEnforceAtomsPeriodic", __FILE__,
            PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
            "Atom %" PetscInt_FMT
            " at (%g, %g, %g), box [%g,%g] x [%g,%g] x [%g,%g], "
            "bcc = (%d,%d,%d)\n",
            site_i, mean_q_i(0), mean_q_i(1), mean_q_i(2), box_coords[0],
            box_coords[3], box_coords[1], box_coords[4], box_coords[2],
            box_coords[5], (int)bcc[0], (int)bcc[1], (int)bcc[2]));
        free(found_local_ptr);
        PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q",
                                      NULL, NULL, (void**)&mean_q_ptr));
        PetscFunctionReturn(PETSC_ERR_RETURN);
      }
    }
  }

  free(found_local_ptr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmSyncCoorFromMeanQ(simulation));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmEnforceGhostAtomsPeriodic(Simulation& simulation) {

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  unsigned int dim = NumberDimensions;

  //! Number of local physical sites
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get new local size
  PetscInt n_sites_local_ghosted;
  PetscCall(
      DMSwarmGetLocalSize(simulation.dm(), &n_sites_local_ghosted));

  //! Get atomistic coordinates
  PetscScalar* mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  //! Get box index field
  PetscInt* box_idx_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "box-idx", NULL, NULL,
                            (void**)&box_idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the global (x,y,z) indices of the lower left corner and size of the
    local region, excluding ghost points.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscReal lmin[3], lmax[3];
  PetscReal gmin[3], gmax[3], mesh_buffer_width;
  DM background_mesh;
  PetscCall(DMSwarmGetCellDM(simulation.dm(), &background_mesh));
  PetscCall(DMGetPhysicalLocalBoundingBox(background_mesh, lmin, lmax));
  PetscCall(DMGetPhysicalBoundingBox(background_mesh, gmin, gmax,
                                     &mesh_buffer_width));

  //! Define local element with buffer width
  const PetscScalar el_coords[6] = {
      lmin[0] - mesh_buffer_width, lmin[1] - mesh_buffer_width,
      lmin[2] - mesh_buffer_width, lmax[0] + mesh_buffer_width,
      lmax[1] + mesh_buffer_width, lmax[2] + mesh_buffer_width};

  const Eigen::Vector3d lattice_x_B(gmax[0] - gmin[0], 0.0, 0.0);
  const Eigen::Vector3d lattice_y_B(0.0, gmax[1] - gmin[1], 0.0);
  const Eigen::Vector3d lattice_z_B(0.0, 0.0, gmax[2] - gmin[2]);

  //! Get mesh boundary conditions
  DMBoundaryType bcc[3] = {DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DM_BOUNDARY_NONE};
  PetscCall(get_mesh_boundary_condition(bcc, &background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over ghost atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = n_sites_local; site_i < n_sites_local_ghosted; site_i++) {

    PetscInt num_sym_box = 27;
    for (unsigned box_idx = 0; box_idx < num_sym_box; box_idx++) {
      int dx, dy, dz;
      MillerIndexToTuple((Miller_Index)box_idx, dx, dy, dz);

      if (!IsDirectionActive(dx, dy, dz, bcc[0], bcc[1], bcc[2])) continue;

      Eigen::Vector3d mean_q_i;
      mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
          mean_q_ptr[site_i * dim + 2];

      Eigen::Matrix3d periodic_matrix_translation = Eigen::Matrix3d::Zero();
      periodic_matrix_translation(0, 0) = dx;
      periodic_matrix_translation(1, 1) = dy;
      periodic_matrix_translation(2, 2) = dz;

      mean_q_i += periodic_matrix_translation *
                  (lattice_x_B + lattice_y_B + lattice_z_B);

      if (In_Out_Mesh(mean_q_i, el_coords) == true) {
        mean_q_ptr[site_i * dim + 0] = mean_q_i(0);
        mean_q_ptr[site_i * dim + 1] = mean_q_i(1);
        mean_q_ptr[site_i * dim + 2] = mean_q_i(2);
        box_idx_ptr[site_i] = box_idx;
        break;
      }
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "box-idx", NULL,
                                NULL, (void**)&box_idx_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode VecEnforceGhostAtomsPeriodic(Vec mean_q,
                                            const PetscInt* box_idx_ptr,
                                            DM background_mesh) {

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  unsigned int dim = NumberDimensions;

  //! Get local size
  PetscInt n_dof_local;
  PetscCall(VecGetLocalSize(mean_q, &n_dof_local));
  PetscInt n_sites_local = n_dof_local / dim;

  //! Get local ghosted of mean_q
  Vec mean_q_loc;
  PetscCall(VecGhostGetLocalForm(mean_q, &mean_q_loc));

  //! Get local size
  PetscInt n_dof_local_ghosted;
  PetscCall(VecGetLocalSize(mean_q_loc, &n_dof_local_ghosted));
  PetscInt n_sites_local_ghosted = n_dof_local_ghosted / dim;

  //! Get atomistic coordinates
  PetscScalar* mean_q_ptr;
  PetscCall(VecGetArray(mean_q_loc, &mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the global (x,y,z) indices of the lower left corner and size of the
    local region, excluding ghost points.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscReal gmin[3], gmax[3];
  PetscCall(DMGetPhysicalBoundingBox(background_mesh, gmin, gmax, NULL));

  const Eigen::Vector3d lattice_x_B(gmax[0] - gmin[0], 0.0, 0.0);
  const Eigen::Vector3d lattice_y_B(0.0, gmax[1] - gmin[1], 0.0);
  const Eigen::Vector3d lattice_z_B(0.0, 0.0, gmax[2] - gmin[2]);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over ghost atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = n_sites_local; site_i < n_sites_local_ghosted; site_i++) {

    const unsigned box_idx_i = box_idx_ptr[site_i];

    if (box_idx_i > 0) {
      int dx, dy, dz;
      MillerIndexToTuple((Miller_Index)box_idx_i, dx, dy, dz);

      Eigen::Vector3d mean_q_i;
      mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
          mean_q_ptr[site_i * dim + 2];

      Eigen::Matrix3d periodic_matrix_translation = Eigen::Matrix3d::Zero();
      periodic_matrix_translation(0, 0) = dx;
      periodic_matrix_translation(1, 1) = dy;
      periodic_matrix_translation(2, 2) = dz;

      mean_q_i += periodic_matrix_translation *
                  (lattice_x_B + lattice_y_B + lattice_z_B);

      mean_q_ptr[site_i * dim + 0] = mean_q_i(0);
      mean_q_ptr[site_i * dim + 1] = mean_q_i(1);
      mean_q_ptr[site_i * dim + 2] = mean_q_i(2);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArray(mean_q_loc, &mean_q_ptr));
  PetscCall(VecGhostRestoreLocalForm(mean_q, &mean_q_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmFixMeanPositionBox(Simulation& simulation, PetscInt FixLabel,
                                         const PetscScalar box_coords[6]) {

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (FixLabel == 0) {
    PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "DMSwarmFixMeanPositionBox",
                         __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                         "Wrong FixLabel value"));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  unsigned int dim = NumberDimensions;

  //! Number of local physical sites
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get atomistic coordinates
  PetscScalar* mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  //! Get fix-mean-q index field
  PetscInt* idx_bcc_mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx-bcc-mean-q", NULL,
                            NULL, (void**)&idx_bcc_mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = 0; site_i < n_sites_local; site_i++) {

    //! Get site position
    Eigen::Vector3d mean_q_i;
    mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
        mean_q_ptr[site_i * dim + 2];

    // Check if inside domain
    if (In_Out_Mesh_Closed(mean_q_i, box_coords) == true) {
      idx_bcc_mean_q_ptr[site_i] = FixLabel;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-bcc-mean-q",
                                NULL, NULL, (void**)&idx_bcc_mean_q_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmApplyDisplacement(Simulation& simulation, PetscInt FixLabel,
                                        PetscScalar displacement_x,
                                        PetscScalar displacement_y,
                                        PetscScalar displacement_z) {

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Get auxiliar data.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (FixLabel == 0) {
    PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "DMSwarmApplyDisplacement",
                         __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                         "Wrong FixLabel value"));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get system topology
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //! Get local number of sites in the simulation (without ghost)
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get local number of sites in the simulation (with ghost)
  PetscInt n_sites_local_ghosted;
  PetscCall(
      DMSwarmGetLocalSize(simulation.dm(), &n_sites_local_ghosted));

  //! Get number of ghost particles
  PetscInt n_sites_ghost = n_sites_local_ghosted - n_sites_local;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Periodic box index
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt* box_idx_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "box-idx", NULL, NULL,
                            (void**)&box_idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index of the particles
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt* idx_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx", NULL, NULL,
                            (void**)&idx_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index for the ghost particles
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof_local = dim * n_sites_local;
  PetscInt n_dof_ghost = dim * n_sites_ghost;
  PetscInt* idx_dof_ghost = (PetscInt*)malloc(n_dof_ghost * sizeof(PetscInt));
  for (int i = 0; i < n_sites_ghost; i++) {
    for (int j = 0; j < dim; j++) {
      idx_dof_ghost[i * dim + j] = idx_q_ptr[n_sites_local + i] * dim + j;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get fix-mean-q index field
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt* idx_bcc_mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx-bcc-mean-q", NULL,
                            NULL, (void**)&idx_bcc_mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get mean position pointer
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar* mean_q_ptr;  //! Mean position pointer
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Loop over atoms
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (int site_i = 0; site_i < n_sites_local; site_i++) {

    //! Get site index
    PetscInt idx_bcc_mean_q_i = idx_bcc_mean_q_ptr[site_i];

    // Check idex
    if (idx_bcc_mean_q_i == FixLabel) {
      // Update site position
      mean_q_ptr[site_i * dim + 0] += displacement_x;
      mean_q_ptr[site_i * dim + 1] += displacement_y;
      mean_q_ptr[site_i * dim + 2] += displacement_z;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get finite element mesh
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DM background_mesh;
  PetscCall(DMSwarmGetCellDM(simulation.dm(), &background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get mean position vector update ghost atoms and enforce periodic bcc
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec X_mean_q;
  PetscCall(VecCreateGhostWithArray(PETSC_COMM_WORLD, n_dof_local,
                                    PETSC_DETERMINE, n_dof_ghost, idx_dof_ghost,
                                    mean_q_ptr, &X_mean_q));
  PetscCall(VecSetUp(X_mean_q));
  PetscCall(VecGhostUpdateBegin(X_mean_q, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(X_mean_q, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(
      VecEnforceGhostAtomsPeriodic(X_mean_q, box_idx_ptr, background_mesh));
  PetscCall(VecDestroy(&X_mean_q));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore mean-q data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmSyncCoorFromMeanQ(simulation));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore Periodic box index
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "box-idx", NULL,
                                NULL, (void**)&box_idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore idx data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx", NULL, NULL,
                                (void**)&idx_q_ptr));
  free(idx_dof_ghost);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Restore fix-mean-q index field
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-bcc-mean-q",
                                NULL, NULL, (void**)&idx_bcc_mean_q_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode VecFixMeanPositionRHS(Vec RHS, Vec mean_q, Vec mean_q_ref,
                                     const PetscInt* idx_bcc_mean_q) {

  PetscFunctionBeginUser;

  //! @brief Auxiliar atomistic variables
  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get local size
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof;
  PetscCall(VecGetLocalSize(mean_q, &n_dof));
  PetscInt n_sites = n_dof / dim;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to the local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec RHS_loc, mean_q_loc, mean_q_ref_loc;
  PetscCall(VecGhostGetLocalForm(RHS, &RHS_loc));
  PetscCall(VecGhostGetLocalForm(mean_q, &mean_q_loc));
  PetscCall(VecGhostGetLocalForm(mean_q_ref, &mean_q_ref_loc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to raw data from PETSc vec
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar* RHS_ptr;
  PetscCall(VecGetArray(RHS_loc, &RHS_ptr));

  const PetscScalar* mean_q_ptr;
  PetscCall(VecGetArrayRead(mean_q_loc, &mean_q_ptr));

  const PetscScalar* mean_q_ref_ptr;
  PetscCall(VecGetArrayRead(mean_q_ref_loc, &mean_q_ref_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = 0; site_i < n_sites; site_i++) {
    if (idx_bcc_mean_q[site_i] != 0) {
      for (PetscInt alpha = 0; alpha < dim; alpha++) {
        RHS_ptr[site_i * dim + alpha] = mean_q_ptr[site_i * dim + alpha] -
                                        mean_q_ref_ptr[site_i * dim + alpha];
      }
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore raw data to PETSc vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArray(RHS_loc, &RHS_ptr));
  PetscCall(VecRestoreArrayRead(mean_q_loc, &mean_q_ptr));
  PetscCall(VecRestoreArrayRead(mean_q_ref_loc, &mean_q_ref_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore local version of the vectors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostRestoreLocalForm(RHS, &RHS_loc));
  PetscCall(VecGhostRestoreLocalForm(mean_q, &mean_q_loc));
  PetscCall(VecGhostRestoreLocalForm(mean_q_ref, &mean_q_ref_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmFixStdvPositionBox(Simulation& simulation, PetscInt FixLabel,
                                         const PetscScalar box_coords[6]) {

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (FixLabel == 0) {
    PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "DMSwarmFixStdvPositionBox",
                         __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                         "Wrong FixLabel value"));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  unsigned int dim = NumberDimensions;

  //! Number of local physical sites
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get atomistic coordinates
  PetscScalar* mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  //! Get fix-stdv-q index field
  PetscInt* idx_bcc_stdv_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx-bcc-stdv-q", NULL,
                            NULL, (void**)&idx_bcc_stdv_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = 0; site_i < n_sites_local; site_i++) {

    //! Get site position
    Eigen::Vector3d mean_q_i;
    mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
        mean_q_ptr[site_i * dim + 2];

    // Check if inside domain
    if (In_Out_Mesh_Closed(mean_q_i, box_coords) == true) {
      idx_bcc_stdv_q_ptr[site_i] = FixLabel;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-bcc-stdv-q",
                                NULL, NULL, (void**)&idx_bcc_stdv_q_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode VecFixStdvPositionRHS(Vec RHS, Vec stdv_q, Vec stdv_q_ref,
                                     const PetscInt* idx_bcc_stdv_q) {

  PetscFunctionBeginUser;

  //! @brief Auxiliar atomistic variables
  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get local size
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_sites;
  PetscCall(VecGetLocalSize(stdv_q, &n_sites));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to the local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec RHS_loc, stdv_q_loc, stdv_q_ref_loc;
  PetscCall(VecGhostGetLocalForm(RHS, &RHS_loc));
  PetscCall(VecGhostGetLocalForm(stdv_q, &stdv_q_loc));
  PetscCall(VecGhostGetLocalForm(stdv_q_ref, &stdv_q_ref_loc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to raw data from PETSc vec
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar* RHS_ptr;
  PetscCall(VecGetArray(RHS_loc, &RHS_ptr));

  const PetscScalar* stdv_q_ptr;
  PetscCall(VecGetArrayRead(stdv_q_loc, &stdv_q_ptr));

  const PetscScalar* stdv_q_ref_ptr;
  PetscCall(VecGetArrayRead(stdv_q_ref_loc, &stdv_q_ref_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = 0; site_i < n_sites; site_i++) {
    if (idx_bcc_stdv_q[site_i] != 0) {
      RHS_ptr[site_i] = stdv_q_ptr[site_i] - stdv_q_ref_ptr[site_i];
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore raw data to PETSc vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArray(RHS_loc, &RHS_ptr));
  PetscCall(VecRestoreArrayRead(stdv_q_loc, &stdv_q_ptr));
  PetscCall(VecRestoreArrayRead(stdv_q_ref_loc, &stdv_q_ref_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore local version of the vectors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostRestoreLocalForm(RHS, &RHS_loc));
  PetscCall(VecGhostRestoreLocalForm(stdv_q, &stdv_q_loc));
  PetscCall(VecGhostRestoreLocalForm(stdv_q_ref, &stdv_q_ref_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmFixChemicalMultiplierBox(
    Simulation& simulation, const PetscScalar box_coords[6]) {

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  unsigned int dim = NumberDimensions;

  //! Number of local physical sites
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get atomistic coordinates
  PetscScalar* mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  //! Get fix chemical multiplier index field
  PetscInt* idx_bcc_gamma_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx-bcc-gamma", NULL,
                            NULL, (void**)&idx_bcc_gamma_ptr));

  //!  Get the indx of the diffusive atoms
  PetscInt* idx_diff_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx-diff", NULL, NULL,
                            (void**)&idx_diff_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = 0; site_i < n_sites_local; site_i++) {
    if (idx_diff_ptr[site_i] == 1) {
      //! Get site position
      Eigen::Vector3d mean_q_i;
      mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
          mean_q_ptr[site_i * dim + 2];

      // Check if inside domain
      if (In_Out_Mesh_Closed(mean_q_i, box_coords) == true) {
        idx_bcc_gamma_ptr[site_i] = 1;
      }
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-bcc-gamma",
                                NULL, NULL, (void**)&idx_bcc_gamma_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore the indx of the diffusive atoms
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-diff", NULL,
                                NULL, (void**)&idx_diff_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmFixThermalMultiplierBox(Simulation& simulation,
                                              const PetscScalar box_coords[6]) {

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  unsigned int dim = NumberDimensions;

  //! Number of local physical sites
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get atomistic coordinates
  PetscScalar* mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  //! Get fix thermal multiplier index field
  PetscInt* idx_bcc_beta_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx-bcc-beta", NULL,
                            NULL, (void**)&idx_bcc_beta_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Loop over atoms
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_i = 0; site_i < n_sites_local; site_i++) {

    //! Get site position
    Eigen::Vector3d mean_q_i;
    mean_q_i << mean_q_ptr[site_i * dim + 0], mean_q_ptr[site_i * dim + 1],
        mean_q_ptr[site_i * dim + 2];

    // Check if inside domain
    if (In_Out_Mesh_Closed(mean_q_i, box_coords) == true) {
      idx_bcc_beta_ptr[site_i] = 1;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore auxiliar data.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-bcc-beta",
                                NULL, NULL, (void**)&idx_bcc_beta_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

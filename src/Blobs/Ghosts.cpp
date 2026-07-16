/**
 * @file Atoms/Ghosts.cpp
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2024-08-22
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "Macros.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "Mesh/Coordinates.hpp"
#include "Mesh/InOut-Mesh.hpp"
#include "petscdmswarm.h"
#include <map>
#include <petscsystypes.h>
#include <vector>

extern PetscMPIInt size_MPI;
extern PetscMPIInt rank_MPI;

typedef struct {
  PetscInt n_ghost_sites;
  PetscInt n_duplicate_groups;
  PetscInt n_redundant_sites;
} GhostIdxUniquenessStats;

static PetscErrorCode
DMSwarmCheckGhostIdxUniqueness(Simulation &simulation,
                               GhostIdxUniquenessStats *stats);

/********************************************************************************/

PetscErrorCode DMSwarmCreateGhostBlobs(Simulation& simulation, double buffer_width) {

  unsigned int dim = NumberDimensions;

  PetscFunctionBeginUser;

  PetscInt* idx_ptr;
  PetscInt* swarm_rank_ptr;
  PetscInt* ghost_ptr;
  PetscInt* rank_ptr;
  PetscInt* box_idx_ptr;
  PetscScalar* mean_q_ptr;
  PetscScalar* mean_p_ptr;
  PetscScalar* rho_ptr;
  PetscScalar* beta_ptr;
  PetscScalar* mass_ptr;
  PetscInt* idx_bcc_mean_q_ptr;
  PetscInt* idx_beta_bcc_ptr;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Check we have the correct migration algorithm
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMSwarmMigrateType atoms_migrate_type;
  DMSwarmGetMigrateType(simulation.dm(), &atoms_migrate_type);

  if (atoms_migrate_type != DMSWARM_MIGRATE_BASIC) {

    PetscCall(
        PetscError(PETSC_COMM_WORLD, __LINE__, "DMSwarmCreateGhostBlobs",
                   __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                   "Set DMSWARM_MIGRATE_BASIC using DMSwarmSetMigrateType"));

    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  //! Get global size
  PetscInt n_global_size = 0;
  PetscCall(DMSwarmGetSize(simulation.dm(), &n_global_size));

  //! Get local size
  PetscInt n_local_size = 0;
  PetscCall(DMSwarmGetLocalSize(simulation.dm(), &n_local_size));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the global (x,y,z) indices of the lower left corner and size of the
    local region, excluding ghost points.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscReal lmin[3], lmax[3];
  DM background_mesh;
  PetscCall(DMSwarmGetCellDM(simulation.dm(), &background_mesh));
  PetscCall(DMGetPhysicalLocalBoundingBox(background_mesh, lmin, lmax));

  //! Lower left coordinates of the brick
  PetscReal X_ll, Y_ll, Z_ll;
  X_ll = lmin[0];
  Y_ll = lmin[1];
  Z_ll = lmin[2];

  //! Upper right coordinates of the brick
  PetscReal X_ur, Y_ur, Z_ur;
  X_ur = lmax[0];
  Y_ur = lmax[1];
  Z_ur = lmax[2];

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Get mesh boundary conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMBoundaryType bcc[3] = {DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DM_BOUNDARY_NONE};
  PetscCall(get_mesh_boundary_condition(bcc, &background_mesh));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the list of processes neighboring this one.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const PetscMPIInt* neig_ranks_ptr;
  PetscInt nranks;
  PetscCall(DMGetNeighbors(background_mesh, &nranks, &neig_ranks_ptr));
  Eigen::Map<const List1D> ranks_list_1D(neig_ranks_ptr, nranks);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Mark ghost particles in the buffer region of the domain.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));
  Eigen::Map<MatrixType> mean_q(mean_q_ptr, n_local_size, dim);

  int num_ghost = 0;
  int* idx_ghost = (int*)malloc(n_global_size * sizeof(int));
  int* rank_ghost = (int*)malloc(n_global_size * sizeof(int));

  for (int site_i = 0; site_i < n_local_size; site_i++) {

    Eigen::Vector3d mean_q_i = mean_q.row(site_i);

    // Direction (-1 -1 -1): n0
    int rank_m1_m1_m1 = ranks_list_1D(0);
    if (rank_m1_m1_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_m1_m1;
        num_ghost++;
      }
    }

    // Direction (0 -1 -1): n1
    int rank_0_m1_m1 = ranks_list_1D(1);
    if (rank_0_m1_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_m1_m1;
        num_ghost++;
      }
    }

    // Direction (1 -1 -1): n2
    int rank_p1_m1_m1 = ranks_list_1D(2);
    if (rank_p1_m1_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_m1_m1;
        num_ghost++;
      }
    }

    // Direction (-1 0 -1): n3
    int rank_m1_0_m1 = ranks_list_1D(3);
    if (rank_m1_0_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_0_m1;
        num_ghost++;
      }
    }

    // Direction (0 0 -1): n4
    int rank_0_0_m1 = ranks_list_1D(4);
    if (rank_0_0_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_0_m1;
        num_ghost++;
      }
    }

    // Direction (1 0 -1): n5
    int rank_p1_0_m1 = ranks_list_1D(5);
    if (rank_p1_0_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_0_m1;
        num_ghost++;
      }
    }

    // Direction (-1 1 -1): n6
    int rank_m1_p1_m1 = ranks_list_1D(6);
    if (rank_m1_p1_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_p1_m1;
        num_ghost++;
      }
    }

    // Direction (0 1 -1): n7
    int rank_0_p1_m1 = ranks_list_1D(7);
    if (rank_0_p1_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_p1_m1;
        num_ghost++;
      }
    }

    // Direction (1 1 -1): n8
    int rank_p1_p1_m1 = ranks_list_1D(8);
    if (rank_p1_p1_m1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ll + buffer_width;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_p1_m1;
        num_ghost++;
      }
    }

    // Direction (-1 -1 0): n9
    int rank_m1_m1_0 = ranks_list_1D(9);
    if (rank_m1_m1_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_m1_0;
        num_ghost++;
      }
    }

    // Direction (0 -1 0): n10
    int rank_0_m1_0 = ranks_list_1D(10);
    if (rank_0_m1_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_m1_0;
        num_ghost++;
      }
    }

    // Direction (1 -1 0): n11
    int rank_p1_m1_0 = ranks_list_1D(11);
    if (rank_p1_m1_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_m1_0;
        num_ghost++;
      }
    }

    // Direction (-1 0 0): n12
    int rank_m1_0_0 = ranks_list_1D(12);
    if (rank_m1_0_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_0_0;
        num_ghost++;
      }
    }

    // Direction (0 0 0) <- current rank (we skip it): n13

    // Direction (1 0 0): n14
    int rank_p1_0_0 = ranks_list_1D(14);
    if (rank_p1_0_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_0_0;
        num_ghost++;
      }
    }

    // Direction (-1 1 0): n15
    int rank_m1_p1_0 = ranks_list_1D(15);
    if (rank_m1_p1_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_p1_0;
        num_ghost++;
      }
    }

    // Direction (0 1 0): n16
    int rank_0_p1_0 = ranks_list_1D(16);
    if (rank_0_p1_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_p1_0;
        num_ghost++;
      }
    }

    // Direction (1 1 0): n17
    int rank_p1_p1_0 = ranks_list_1D(17);
    if (rank_p1_p1_0 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ll;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_p1_0;
        num_ghost++;
      }
    }

    // Direction (-1 -1 1): n18
    int rank_m1_m1_p1 = ranks_list_1D(18);
    if (rank_m1_m1_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_m1_p1;
        num_ghost++;
      }
    }

    // Direction (0 -1 1): n19
    int rank_0_m1_p1 = ranks_list_1D(19);
    if (rank_0_m1_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_m1_p1;
        num_ghost++;
      }
    }

    // Direction (1 -1 1): n20
    int rank_p1_m1_p1 = ranks_list_1D(20);
    if (rank_p1_m1_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ll + buffer_width;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_m1_p1;
        num_ghost++;
      }
    }

    // Direction (-1 0 1): n21
    int rank_m1_0_p1 = ranks_list_1D(21);
    if (rank_m1_0_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_0_p1;
        num_ghost++;
      }
    }

    // Direction (0 0 1): n22
    int rank_0_0_p1 = ranks_list_1D(22);
    if (rank_0_0_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_0_p1;
        num_ghost++;
      }
    }

    // Direction (1 0 1): n23
    int rank_p1_0_p1 = ranks_list_1D(23);
    if (rank_p1_0_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ll;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_0_p1;
        num_ghost++;
      }
    }

    // Direction (-1 1 1): n24
    int rank_m1_p1_p1 = ranks_list_1D(24);
    if (rank_m1_p1_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ll + buffer_width;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_m1_p1_p1;
        num_ghost++;
      }
    }

    // Direction (0 1 1): n25
    int rank_0_p1_p1 = ranks_list_1D(25);
    if (rank_0_p1_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ll;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_0_p1_p1;
        num_ghost++;
      }
    }

    // Direction (1 1 1): n26
    int rank_p1_p1_p1 = ranks_list_1D(26);
    if (rank_p1_p1_p1 >= 0) {

      PetscScalar buffer_coords[6];
      buffer_coords[0] = X_ur - buffer_width;
      buffer_coords[1] = Y_ur - buffer_width;
      buffer_coords[2] = Z_ur - buffer_width;

      buffer_coords[3] = X_ur;
      buffer_coords[4] = Y_ur;
      buffer_coords[5] = Z_ur;

      if (In_Out_Mesh_Closed(mean_q_i, buffer_coords)) {  //!
        idx_ghost[num_ghost] = site_i;
        rank_ghost[num_ghost] = rank_p1_p1_p1;
        num_ghost++;
      }
    }
  }

  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL, NULL,
                                (void**)&mean_q_ptr));

#ifdef DEBUG_MODE
  std::cout << "Number of ghost atoms: " << num_ghost
            << " from rank: " << rank_MPI << std::endl;
#endif

  PetscCall(DMSwarmAddNPoints(simulation.dm(), num_ghost));

  //! Initialize new points
  PetscCall(DMSwarmGetField(simulation.dm(), "idx", NULL, NULL,
                            (void**)&idx_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), DMSwarmField_rank, NULL,
                            NULL, (void**)&swarm_rank_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "ghost", NULL, NULL,
                            (void**)&ghost_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "box-idx", NULL, NULL,
                            (void**)&box_idx_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "MPI-rank", NULL, NULL,
                            (void**)&rank_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "mean-p", NULL, NULL,
                            (void**)&mean_p_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "rho", NULL, NULL,
                            (void**)&rho_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "beta", NULL, NULL,
                            (void**)&beta_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "mass", NULL, NULL,
                            (void**)&mass_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "idx-bcc-mean-q", NULL,
                            NULL, (void**)&idx_bcc_mean_q_ptr));

  PetscCall(DMSwarmGetField(simulation.dm(), "idx-bcc-beta", NULL,
                            NULL, (void**)&idx_beta_bcc_ptr));

  for (int idx = 0; idx < num_ghost; idx++) {

    int loc_site_i = idx_ghost[idx];
    int ghost_i = n_local_size + idx;

    swarm_rank_ptr[ghost_i] = rank_ghost[idx];
    ghost_ptr[ghost_i] = 1;
    box_idx_ptr[ghost_i] = box_idx_ptr[loc_site_i];
    rank_ptr[ghost_i] = rank_MPI;
    idx_ptr[ghost_i] = idx_ptr[loc_site_i];
    for (int alpha = 0; alpha < dim; alpha++) {
      mean_q_ptr[ghost_i * dim + alpha] = mean_q_ptr[loc_site_i * dim + alpha];
      mean_p_ptr[ghost_i * dim + alpha] =
          mean_p_ptr[loc_site_i * dim + alpha];
    }
    rho_ptr[ghost_i] = rho_ptr[loc_site_i];
    beta_ptr[ghost_i] = beta_ptr[loc_site_i];
    mass_ptr[ghost_i] = mass_ptr[loc_site_i];
    idx_bcc_mean_q_ptr[ghost_i] = idx_bcc_mean_q_ptr[loc_site_i];
    idx_beta_bcc_ptr[ghost_i] = idx_beta_bcc_ptr[loc_site_i];
  }

  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx", NULL, NULL,
                                (void**)&idx_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL, NULL,
                                (void**)&mean_q_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-p", NULL, NULL,
                                (void**)&mean_p_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "ghost", NULL, NULL,
                                (void**)&ghost_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "box-idx", NULL,
                                NULL, (void**)&box_idx_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "MPI-rank", NULL,
                                NULL, (void**)&rank_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), DMSwarmField_rank,
                                NULL, NULL, (void**)&swarm_rank_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "rho", NULL, NULL,
                                (void**)&rho_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "beta", NULL, NULL,
                                (void**)&beta_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "mass", NULL, NULL,
                                (void**)&mass_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-bcc-mean-q",
                                NULL, NULL, (void**)&idx_bcc_mean_q_ptr));

  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx-bcc-beta",
                                NULL, NULL, (void**)&idx_beta_bcc_ptr));

  PetscCall(DMSwarmSyncCoorFromMeanQ(simulation));

  //! Migrate new points
  PetscCall(DMSwarmMigrate(simulation.dm(), PETSC_TRUE));

  //! Enforce periodic boundary conditions over the ghost blobs
  PetscCall(DMSwarmEnforceGhostBlobsPeriodic(simulation));

  //! Check the uniqueness of the ghost indices
  PetscCall(DMSwarmCheckGhostIdxUniqueness(simulation, NULL));

  //! Destroy auxiliar variables
  free(idx_ghost);
  free(rank_ghost);

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

namespace {

struct LocalSiteRecord {
  PetscInt local_site = 0;
  PetscInt idx = 0;
  PetscInt box_idx = 0;
  PetscInt is_ghost = 0;
};

}  // namespace

/********************************************************************************/

static PetscErrorCode
DMSwarmCheckGhostIdxUniqueness(Simulation &simulation,
                               GhostIdxUniquenessStats *stats) {

  PetscFunctionBeginUser;
  PetscValidHeaderSpecific(simulation.dm(), DM_CLASSID, 1);

  PetscInt *idx_ptr = NULL;
  PetscInt *ghost_ptr = NULL;
  PetscInt *box_idx_ptr = NULL;
  PetscCall(DMSwarmGetField(simulation.dm(), "idx", NULL, NULL,
                            (void **)&idx_ptr));
  PetscCall(DMSwarmGetField(simulation.dm(), "ghost", NULL, NULL,
                            (void **)&ghost_ptr));
  PetscCall(DMSwarmGetField(simulation.dm(), "box-idx", NULL, NULL,
                            (void **)&box_idx_ptr));

  PetscInt n_sites_local_ghosted = 0;
  PetscCall(DMSwarmGetLocalSize(simulation.dm(), &n_sites_local_ghosted));

  std::map<PetscInt, std::vector<LocalSiteRecord>> sites_by_idx;

  for (PetscInt site_i = 0; site_i < n_sites_local_ghosted; site_i++) {
    LocalSiteRecord record;
    record.local_site = site_i;
    record.idx = idx_ptr[site_i];
    record.box_idx = box_idx_ptr[site_i];
    record.is_ghost = ghost_ptr[site_i];
    sites_by_idx[record.idx].push_back(record);
  }

  PetscInt n_duplicate_groups_local = 0;
  PetscInt n_redundant_sites_local = 0;
  PetscInt n_ghost_sites_local = 0;

  for (PetscInt site_i = 0; site_i < n_sites_local_ghosted; site_i++) {
    if (ghost_ptr[site_i] != 0) {
      n_ghost_sites_local++;
    }
  }

  for (const auto &entry : sites_by_idx) {
    const PetscInt idx = entry.first;
    const std::vector<LocalSiteRecord> &records = entry.second;
    if (records.size() <= 1) {
      continue;
    }

    n_duplicate_groups_local++;
    n_redundant_sites_local += static_cast<PetscInt>(records.size()) - 1;

    for (size_t k = 1; k < records.size(); k++) {
      PetscCall(PetscPrintf(
          PETSC_COMM_WORLD,
          "[Ghost-MPI-uniqueness] duplicate swarm idx %" PetscInt_FMT
          "\n"
          "  reference local %" PetscInt_FMT " ghost=%" PetscInt_FMT
          " box-idx %" PetscInt_FMT "\n"
          "  duplicate local %" PetscInt_FMT " ghost=%" PetscInt_FMT
          " box-idx %" PetscInt_FMT "\n",
          idx, records.front().local_site, records.front().is_ghost,
          records.front().box_idx, records[k].local_site, records[k].is_ghost,
          records[k].box_idx));
    }
  }

  PetscInt n_duplicate_groups_global = 0;
  PetscInt n_redundant_sites_global = 0;
  PetscInt n_ghost_sites_global = 0;

  PetscCall(MPIU_Allreduce(&n_duplicate_groups_local, &n_duplicate_groups_global,
                           1, MPIU_INT, MPIU_SUM, PETSC_COMM_WORLD));
  PetscCall(MPIU_Allreduce(&n_redundant_sites_local, &n_redundant_sites_global,
                           1, MPIU_INT, MPIU_SUM, PETSC_COMM_WORLD));
  PetscCall(MPIU_Allreduce(&n_ghost_sites_local, &n_ghost_sites_global, 1,
                           MPIU_INT, MPIU_SUM, PETSC_COMM_WORLD));

  if (stats != NULL) {
    stats->n_ghost_sites = n_ghost_sites_global;
    stats->n_duplicate_groups = n_duplicate_groups_global;
    stats->n_redundant_sites = n_redundant_sites_global;
  }

  PetscCall(DMSwarmRestoreField(simulation.dm(), "box-idx", NULL, NULL,
                                (void **)&box_idx_ptr));
  PetscCall(DMSwarmRestoreField(simulation.dm(), "ghost", NULL, NULL,
                                (void **)&ghost_ptr));
  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx", NULL, NULL,
                                (void **)&idx_ptr));

  PetscCheck(n_duplicate_groups_global == 0, PETSC_COMM_WORLD, PETSC_ERR_USER,
             "Ghost MPI idx uniqueness check failed: %" PetscInt_FMT
             " duplicate idx groups",
             n_duplicate_groups_global);

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmDestroyGhostBlobs(Simulation& simulation) {

  PetscFunctionBeginUser;

  //! Get global size
  PetscInt n_global_size = 0;
  PetscCall(DMSwarmGetSize(simulation.dm(), &n_global_size));

  //! Get local size
  PetscInt n_local_size = 0;
  PetscCall(DMSwarmGetLocalSize(simulation.dm(), &n_local_size));

  //!
  PetscInt num_remove_point = 0;

  PetscInt* swarm_rank_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), DMSwarmField_rank, NULL,
                            NULL, (void**)&swarm_rank_ptr));

  PetscInt* rank_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "MPI-rank", NULL, NULL,
                            (void**)&rank_ptr));

  for (int site_i = 0; site_i < n_local_size; site_i++) {

    if ((swarm_rank_ptr[site_i] != rank_ptr[site_i]) &&
        (swarm_rank_ptr[site_i] == rank_MPI)) {
      num_remove_point += 1;
    }

    if ((swarm_rank_ptr[site_i] != rank_ptr[site_i]) &&
        (rank_ptr[site_i] != rank_MPI)) {
      swarm_rank_ptr[site_i] = rank_ptr[site_i];
    }
  }

  PetscCall(DMSwarmRestoreField(simulation.dm(), DMSwarmField_rank,
                                NULL, NULL, (void**)&swarm_rank_ptr));
  PetscCall(DMSwarmRestoreField(simulation.dm(), "MPI-rank", NULL,
                                NULL, (void**)&rank_ptr));

  for (int site_i = 0; site_i < num_remove_point; site_i++) {
    PetscCall(DMSwarmRemovePoint(simulation.dm()));
  }

#ifdef DEBUG_MODE
  std::cout << "Number of removed ghost atoms: " << num_remove_point
            << " at rank: " << rank_MPI << std::endl;
#endif

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

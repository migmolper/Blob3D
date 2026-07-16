/**
 * @file Blobs/Neighbors.cpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief Create and destroy the neighbors list for the particles
 */

#include <cstdlib>
#if __APPLE__
#include <malloc/_malloc.h>
#endif
#include <petscsystypes.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "Blobs/Blobs.hpp"
#include "Blobs/Neighbors.hpp"
#include "Macros.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "petscdmswarm.h"
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
extern PetscMPIInt size_MPI;
extern PetscMPIInt rank_MPI;

/************************************************************************/

PetscErrorCode DMSwarmCreateNeighborsBlobs(Simulation &simulation,
                                           double r_cutoff) {

  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  //! Get the local size
  PetscInt n_sites_local = simulation.n_sites_local();

  //! Get the local size (ghosted)
  PetscInt n_sites_local_ghosted;
  PetscCall(DMSwarmGetLocalSize(simulation.dm(), &n_sites_local_ghosted));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get mean position vector
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *mean_q_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "mean-q", NULL, NULL,
                            (void **)&mean_q_ptr));
  Eigen::Map<MatrixType> mean_q(mean_q_ptr, n_sites_local_ghosted, 3);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index of the particles
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt *idx_ptr;
  PetscCall(
      DMSwarmGetField(simulation.dm(), "idx", NULL, NULL, (void **)&idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get ghost indicator
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt *ghost_ptr;
  PetscCall(DMSwarmGetField(simulation.dm(), "ghost", NULL, NULL,
                            (void **)&ghost_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Find neighbors
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  IS *mechanical_neighs_idx;
  PetscCall(PetscMalloc1(n_sites_local_ghosted, &mechanical_neighs_idx));

  for (PetscInt site_i = 0; site_i < n_sites_local_ghosted; site_i++) {

    //! Get mean position of site i
    Eigen::Vector3d mean_q_i = mean_q.block<1, 3>(site_i, 0);

    //! Create auxiliar list to store the diffusive neighbors of the site i
    unsigned int num_neigh_i = 0;
    PetscInt *mech_neighs_idx_i_ptr =
        (PetscInt *)calloc(maxneigh, sizeof(PetscInt));

    //! Add in the first position the site i
    mech_neighs_idx_i_ptr[0] = site_i;
    num_neigh_i += 1;

    //! Search neighbourhs
    for (unsigned site_j = 0; site_j < n_sites_local_ghosted; site_j++) {

      if (idx_ptr[site_i] != idx_ptr[site_j]) {

        //! Get mean position and specie of site j in the periodic box
        Eigen::Vector3d mean_q_j = mean_q.block<1, 3>(site_j, 0);

        //! Check if site j is the neibourhood of the site i
        double norm_r_ij = (mean_q_i - mean_q_j).norm();

        //!
        if ((norm_r_ij <= r_cutoff) && (num_neigh_i < maxneigh)) {

          mech_neighs_idx_i_ptr[num_neigh_i] = site_j;
          num_neigh_i += 1;

        } else if ((norm_r_ij <= r_cutoff) && (num_neigh_i == maxneigh)) {

          PetscCall(PetscError(
              PETSC_COMM_SELF, __LINE__, "DMSwarmCreateNeighborsBlobs",
              __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
              "Max number of neighbors reached for particle %i", site_i));
          PetscFunctionReturn(PETSC_ERR_RETURN);
        }
      }
    }

    if ((ghost_ptr[site_i] == 0) && (num_neigh_i == 1)) {
      PetscCall(PetscPrintf(PETSC_COMM_SELF,
                            "No neighbors where finded for particle %i at "
                            "(%f, %f, %f). Increase the search radius.\n",
                            idx_ptr[site_i], mean_q_i(0), mean_q_i(1),
                            mean_q_i(2)));
    }

    //! Create list of neighbors for the site i
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, num_neigh_i,
                              mech_neighs_idx_i_ptr, PETSC_COPY_VALUES,
                              &mechanical_neighs_idx[site_i]));

    free(mech_neighs_idx_i_ptr);
  }

  //! Wait
  PetscCall(PetscBarrier(NULL));

  //! Transfer ownership of neighbor lists to the simulation topology
  simulation.topology().adopt(mechanical_neighs_idx, n_sites_local_ghosted);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore data
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "mean-q", NULL, NULL,
                                (void **)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore idx data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "idx", NULL, NULL,
                                (void **)&idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore ghost indicator
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(simulation.dm(), "ghost", NULL, NULL,
                                (void **)&ghost_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmDestroyNeighborsBlobs(Simulation &simulation) {
  PetscFunctionBeginUser;
  PetscCall(simulation.topology().clear());
  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmGetParticleNeighbors(ParticleTopology *particle_topology,
                                           IS mechanical_neighs_idx) {

  PetscFunctionBeginUser;

  PetscCall(ISGetSize(mechanical_neighs_idx, &particle_topology->size));

  //! Get the indices of the mechanical neighbors
  if (particle_topology->size > 0) {
    PetscCall(ISGetIndices(mechanical_neighs_idx, &particle_topology->list));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode
DMSwarmRestoreParticleNeighbors(ParticleTopology *particle_topology,
                                IS mechanical_neighs_idx) {

  PetscFunctionBeginUser;

  if (particle_topology->size > 0) {
    PetscCall(
        ISRestoreIndices(mechanical_neighs_idx, &particle_topology->list));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/
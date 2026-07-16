/**
 * @file Atoms/Topology.cpp
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2024-08-02
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "Blobs/Ghosts.hpp"
#include "Blobs/Neighbors.hpp"
#include "Macros.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "Mesh/Update-Mesh.hpp"
#include "petscis.h"
#include <petscerror.h>
#include <petscsystypes.h>

extern PetscMPIInt size_MPI;
extern PetscMPIInt rank_MPI;

/********************************************************************************/

PetscErrorCode DMSwarmGenerateBlobsTopology(DMD* Simulation,
                                            double buffer_width) {

  PetscFunctionBeginUser;

  //! 1º Create ghost atoms
  PetscCall(
      DMSwarmSetMigrateType(Simulation->atomistic_data, DMSWARM_MIGRATE_BASIC));
  PetscCall(DMSwarmCreateGhostBlobs(Simulation, buffer_width));

  //! 2º Compute list of neighbors
  PetscCall(DMSwarmCreateNeighborsBlobs(Simulation, buffer_width));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmRegenerateBlobsTopology(DMD* Simulation,
                                              double buffer_width,
                                              PetscBool MIGRATE_BLOBS,
                                              PetscBool PARTICLE_INSERTION) {

  PetscFunctionBeginUser;

  DM FE_Mesh;
  PetscInt n_sites_global = 0;

  //! 1º: Destroy topology
  PetscCall(DMSwarmDestroyNeighborsBlobs(Simulation));

  //! 2º: Destroy ghost blobs
  if (MIGRATE_BLOBS == PETSC_TRUE) {
    PetscCall(DMSwarmDestroyGhostBlobs(Simulation));
  }

  //! 3º: Rebin blobs and update number of blobs and check consistency
  if (MIGRATE_BLOBS == PETSC_TRUE) {
    PetscCall(DMSwarmSetMigrateType(Simulation->atomistic_data,
                                    DMSWARM_MIGRATE_DMCELLNSCATTER));
    PetscCall(DMSwarmMigrate(Simulation->atomistic_data, PETSC_TRUE));
    PetscCall(DMSwarmGetSize(Simulation->atomistic_data, &n_sites_global));
    PetscCall(DMSwarmGetLocalSize(Simulation->atomistic_data,
                                  &(Simulation->n_sites_local)));
    if (n_sites_global != Simulation->n_sites_global) {
      PetscCall(PetscError(
          PETSC_COMM_SELF, __LINE__, "DMSwarmRegenerateBlobsTopology", __FILE__,
          PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
          "Number of blobs is not consistent: %i (new), %i (old) \n",
          n_sites_global, Simulation->n_sites_global));

      PetscFunctionReturn(PETSC_ERR_RETURN);
    }
  }

  //! 4º: Insert new particles if needed
/*   if (PARTICLE_INSERTION == PETSC_TRUE) {

    PetscCall(DMSwarmAddNPoints(Simulation->atomistic_data, num_ghost));

    PetscCall(DMSwarmGetLocalSize(Simulation->atomistic_data,
                                  &(Simulation->n_sites_local)));
  } */

  //! 5º: Update "MPI-rank"
  if (MIGRATE_BLOBS == PETSC_TRUE) {
    PetscInt* rank_ptr;
    PetscCall(DMSwarmGetField(Simulation->atomistic_data, "MPI-rank", NULL,
                              NULL, (void**)&rank_ptr));
#pragma omp parallel for schedule(runtime)
    for (PetscInt site_u = 0; site_u < Simulation->n_sites_local; site_u++) {
      rank_ptr[site_u] = rank_MPI;
    }
    PetscCall(DMSwarmRestoreField(Simulation->atomistic_data, "MPI-rank", NULL,
                                  NULL, (void**)&rank_ptr));
  }

  //! 6º: create new ghost blobs and impose periodicity over them
  if (MIGRATE_BLOBS == PETSC_TRUE) {
    PetscCall(DMSwarmSetMigrateType(Simulation->atomistic_data,
                                    DMSWARM_MIGRATE_BASIC));
    PetscCall(DMSwarmCreateGhostBlobs(Simulation, buffer_width));
  }

  //! 7º: Update topology
  PetscCall(DMSwarmCreateNeighborsBlobs(Simulation, buffer_width));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmDestroyBlobsTopology(DMD* Simulation) {

  PetscFunctionBeginUser;

  //! 1º Destroy list of neighs
  PetscCall(DMSwarmDestroyNeighborsBlobs(Simulation));

  //! 2º Destroy ghost blobs
  PetscCall(DMSwarmDestroyGhostBlobs(Simulation));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

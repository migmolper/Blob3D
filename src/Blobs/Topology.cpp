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
#include "Mesh/Coordinates.hpp"
#include "petscis.h"
#include <petscerror.h>
#include <petscsystypes.h>

extern PetscMPIInt size_MPI;
extern PetscMPIInt rank_MPI;

/********************************************************************************/

PetscErrorCode DMSwarmGenerateBlobsTopology(DMD *Simulation,
                                            double buffer_width) {

  PetscFunctionBeginUser;

  //! 1º Create ghost atoms
  PetscCall(
      DMSwarmSetMigrateType(Simulation->atomistic_data, DMSWARM_MIGRATE_BASIC));
  PetscCall(DMSwarmCreateGhostBlobs(Simulation, buffer_width));

  //! 2º Compute list of mechanical neighs
  PetscCall(DMSwarmCreateNeighborsBlobs(Simulation, buffer_width));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmRegenerateBlobsTopology(DMD *Simulation,
                                              double buffer_width) {

  PetscFunctionBeginUser;

  PetscInt n_sites_global = 0;

  //! 1º: Destroy topology
  PetscCall(DMSwarmDestroyNeighborsBlobs(Simulation));

  //! 2º: Destroy ghost atoms
  PetscCall(DMSwarmDestroyGhostBlobs(Simulation));

  //! 3º: Enforce periodic conditions over physical atoms
  PetscCall(DMSwarmEnforceAtomsPeriodic(Simulation, buffer_width));

  //! 4º: Rebin atoms
  PetscCall(DMSwarmSyncCoorFromMeanQ(Simulation));
  PetscCall(DMSwarmSetMigrateType(Simulation->atomistic_data,
                                  DMSWARM_MIGRATE_DMCELLNSCATTER));
  PetscCall(DMSwarmMigrate(Simulation->atomistic_data, PETSC_TRUE));

  //! 5º: Update number of atoms and check consistency
  PetscCall(DMSwarmGetSize(Simulation->atomistic_data, &n_sites_global));
  PetscCall(DMSwarmGetLocalSize(Simulation->atomistic_data,
                                &(Simulation->n_sites_local)));
  if (n_sites_global != Simulation->n_sites_global) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_RETURN,
            "Number of atoms is not consistent: %" PetscInt_FMT
            " (new), %" PetscInt_FMT " (old)",
            n_sites_global, Simulation->n_sites_global);
  }

  //! 6º: Update "MPI-rank"
  PetscInt *rank_ptr;
  PetscCall(DMSwarmGetField(Simulation->atomistic_data, "MPI-rank", NULL, NULL,
                            (void **)&rank_ptr));
  for (PetscInt site_u = 0; site_u < Simulation->n_sites_local; site_u++) {
    rank_ptr[site_u] = rank_MPI;
  }
  PetscCall(DMSwarmRestoreField(Simulation->atomistic_data, "MPI-rank", NULL,
                                NULL, (void **)&rank_ptr));

  //! 7º: create new ghost atoms and impose periodicity over them
  PetscCall(
      DMSwarmSetMigrateType(Simulation->atomistic_data, DMSWARM_MIGRATE_BASIC));
  PetscCall(DMSwarmCreateGhostBlobs(Simulation, buffer_width));

  //! 8º: Update topology
  PetscCall(DMSwarmCreateNeighborsBlobs(Simulation, buffer_width));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode DMSwarmDestroyBlobsTopology(DMD *Simulation) {

  PetscFunctionBeginUser;

  //! 1º Destroy list of neighs
  PetscCall(DMSwarmDestroyNeighborsBlobs(Simulation));

  //! 2º Destroy ghost atoms
  PetscCall(DMSwarmDestroyGhostBlobs(Simulation));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

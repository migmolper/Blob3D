/**
 * @file Simulation.cpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief RAII helpers for Blob3D Simulation.
 * @version 0.1
 * @date 2026-07-16
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "Simulation.hpp"

#include "Blobs/Ghosts.hpp"
#include "Blobs/Neighbors.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "Mesh/Coordinates.hpp"

extern PetscMPIInt rank_MPI;

/********************************************************************************/

/* ParticleSwarm */

/********************************************************************************/

ParticleSwarm::~ParticleSwarm() { release(); }

ParticleSwarm::ParticleSwarm(ParticleSwarm &&other) noexcept
    : dm_(other.dm_), dump2petsc_(other.dump2petsc_),
      n_global_(other.n_global_), n_local_(other.n_local_),
      n_ghost_(other.n_ghost_) {
  other.dm_ = nullptr;
  other.dump2petsc_ = nullptr;
  other.n_global_ = 0;
  other.n_local_ = 0;
  other.n_ghost_ = 0;
}

ParticleSwarm &ParticleSwarm::operator=(ParticleSwarm &&other) noexcept {
  if (this != &other) {
    release();
    dm_ = other.dm_;
    dump2petsc_ = other.dump2petsc_;
    n_global_ = other.n_global_;
    n_local_ = other.n_local_;
    n_ghost_ = other.n_ghost_;
    other.dm_ = nullptr;
    other.dump2petsc_ = nullptr;
    other.n_global_ = 0;
    other.n_local_ = 0;
    other.n_ghost_ = 0;
  }
  return *this;
}

void ParticleSwarm::release() noexcept {
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Destroy the background mesh attached to the swarm
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (dm_) {
    DM fe_mesh = nullptr;
    if (DMSwarmGetCellDM(dm_, &fe_mesh) == PETSC_SUCCESS && fe_mesh) {
      DMDestroy(&fe_mesh);
    }
    DMDestroy(&dm_);
    dm_ = nullptr;
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Destroy dump → PETSc mapping
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (dump2petsc_) {
    AODestroy(&dump2petsc_);
    dump2petsc_ = nullptr;
  }

  n_global_ = 0;
  n_local_ = 0;
  n_ghost_ = 0;
}

void ParticleSwarm::adopt(DM dm, AO dump2petsc, PetscInt n_global,
                          PetscInt n_local) {
  release();
  dm_ = dm;
  dump2petsc_ = dump2petsc;
  n_global_ = n_global;
  n_local_ = n_local;
  n_ghost_ = 0;
}

/********************************************************************************/

/* NeighborTopology */

/********************************************************************************/

NeighborTopology::~NeighborTopology() {
  //! Best-effort cleanup; PETSc may already be finalized
  if (neighs_) {
    for (PetscInt i = 0; i < n_; ++i) {
      if (neighs_[i]) {
        ISDestroy(&neighs_[i]);
      }
    }
    PetscFree(neighs_);
    neighs_ = nullptr;
    n_ = 0;
  }
}

NeighborTopology::NeighborTopology(NeighborTopology &&other) noexcept
    : neighs_(other.neighs_), n_(other.n_) {
  other.neighs_ = nullptr;
  other.n_ = 0;
}

NeighborTopology &
NeighborTopology::operator=(NeighborTopology &&other) noexcept {
  if (this != &other) {
    clear();
    neighs_ = other.neighs_;
    n_ = other.n_;
    other.neighs_ = nullptr;
    other.n_ = 0;
  }
  return *this;
}

void NeighborTopology::adopt(IS *neighs, PetscInt n) {
  clear();
  neighs_ = neighs;
  n_ = n;
}

PetscErrorCode NeighborTopology::clear() {
  PetscFunctionBeginUser;

  if (!neighs_) {
    n_ = 0;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Destroy mechanical neighs
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (PetscInt i = 0; i < n_; ++i) {
    PetscCall(ISDestroy(&neighs_[i]));
  }
  PetscCall(PetscFree(neighs_));
  neighs_ = nullptr;
  n_ = 0;

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

/* Simulation */

/********************************************************************************/

Simulation::~Simulation() {
  //! Destroy topology while the DM is still alive, then release particles
  if (dm()) {
    (void)destroy_topology();
  }
}

/********************************************************************************/

PetscErrorCode Simulation::generate_topology(double buffer_width) {

  PetscFunctionBeginUser;

  //! 1º Create ghost atoms
  PetscCall(DMSwarmSetMigrateType(dm(), DMSWARM_MIGRATE_BASIC));
  PetscCall(DMSwarmCreateGhostBlobs(*this, buffer_width));

  //! 2º Compute list of mechanical neighs
  PetscCall(DMSwarmCreateNeighborsBlobs(*this, buffer_width));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode Simulation::regenerate_topology(double buffer_width) {

  PetscFunctionBeginUser;

  PetscInt n_sites_global = 0;

  //! 1º: Destroy topology
  PetscCall(DMSwarmDestroyNeighborsBlobs(*this));

  //! 2º: Destroy ghost atoms
  PetscCall(DMSwarmDestroyGhostBlobs(*this));

  //! 3º: Enforce periodic conditions over physical atoms
  PetscCall(DMSwarmEnforceAtomsPeriodic(*this, buffer_width));

  //! 4º: Rebin atoms
  PetscCall(DMSwarmSyncCoorFromMeanQ(*this));
  PetscCall(DMSwarmSetMigrateType(dm(), DMSWARM_MIGRATE_DMCELLNSCATTER));
  PetscCall(DMSwarmMigrate(dm(), PETSC_TRUE));

  //! 5º: Update number of particles and check consistency
  PetscCall(DMSwarmGetSize(dm(), &n_sites_global));
  PetscCall(DMSwarmGetLocalSize(dm(), &n_sites_local()));
  if (n_sites_global != this->n_sites_global()) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_RETURN,
            "Number of particles is not consistent: %" PetscInt_FMT
            " (new), %" PetscInt_FMT " (old)",
            n_sites_global, this->n_sites_global());
  }

  //! 6º: Update "MPI-rank"
  PetscInt *rank_ptr;
  PetscCall(DMSwarmGetField(dm(), "MPI-rank", NULL, NULL, (void **)&rank_ptr));
  for (PetscInt site_u = 0; site_u < n_sites_local(); site_u++) {
    rank_ptr[site_u] = rank_MPI;
  }
  PetscCall(
      DMSwarmRestoreField(dm(), "MPI-rank", NULL, NULL, (void **)&rank_ptr));

  //! 7º: create new ghost atoms and impose periodicity over them
  PetscCall(DMSwarmSetMigrateType(dm(), DMSWARM_MIGRATE_BASIC));
  PetscCall(DMSwarmCreateGhostBlobs(*this, buffer_width));

  //! 8º: Update topology
  PetscCall(DMSwarmCreateNeighborsBlobs(*this, buffer_width));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

PetscErrorCode Simulation::destroy_topology() {

  PetscFunctionBeginUser;

  //! 1º Destroy list of neighs
  PetscCall(DMSwarmDestroyNeighborsBlobs(*this));

  //! 2º Destroy ghost atoms
  PetscCall(DMSwarmDestroyGhostBlobs(*this));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

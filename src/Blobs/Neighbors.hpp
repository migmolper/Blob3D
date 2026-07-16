/**
 * @file Atoms/Neighbors.hpp
 */

#ifndef NEIGHBORS_BLOBS_HPP
#define NEIGHBORS_BLOBS_HPP

#include "Simulation.hpp"
#include <Eigen/Dense>

PetscErrorCode DMSwarmCreateNeighborsBlobs(Simulation& simulation,
                                           double r_cutoff);

PetscErrorCode DMSwarmDestroyNeighborsBlobs(Simulation& simulation);

PetscErrorCode DMSwarmGetParticleNeighbors(ParticleTopology* particle_topology,
                                           IS mechanical_neighs_idx);

PetscErrorCode DMSwarmRestoreParticleNeighbors(
    ParticleTopology* particle_topology, IS mechanical_neighs_idx);

#endif /* NEIGHBORS_BLOBS_HPP */

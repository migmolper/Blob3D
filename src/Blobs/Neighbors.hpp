/**
 * @file Blobs/Neighbors.hpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief Create and destroy the neighbors list for the particles
 */

#ifndef NEIGHBORS_BLOBS_HPP
#define NEIGHBORS_BLOBS_HPP

#include "Simulation.hpp"
#include <Eigen/Dense>

PetscErrorCode DMSwarmCreateNeighborsBlobs(Simulation &simulation);

PetscErrorCode DMSwarmDestroyNeighborsBlobs(Simulation &simulation);

PetscErrorCode DMSwarmGetParticleNeighbors(ParticleTopology *particle_topology,
                                           IS mechanical_neighs_idx);

PetscErrorCode
DMSwarmRestoreParticleNeighbors(ParticleTopology *particle_topology,
                                IS mechanical_neighs_idx);

#endif /* NEIGHBORS_BLOBS_HPP */

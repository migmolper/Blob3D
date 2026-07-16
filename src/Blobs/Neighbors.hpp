/**
 * @file Atoms/Neighbors.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef NEIGHBORS_BLOBS_HPP
#define NEIGHBORS_BLOBS_HPP

#include "Macros.hpp"
#include <Eigen/Dense>

/**
 * @brief Function to compute the list of neighbors for the all simulation
 *
 * @param Simulation
 * @param r_cutoff
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmCreateNeighborsBlobs(DMD *Simulation, double r_cutoff);

/**
 * @brief
 *
 * @param Simulation
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmDestroyNeighborsBlobs(DMD *Simulation);

/**
 * @brief
 *
 * @param particle_topology
 * @param mechanical_neighs_idx
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmGetParticleNeighbors(ParticleTopology *particle_topology,
                                           IS mechanical_neighs_idx);

/**
 * @brief
 *
 * @param particle_topology
 * @param mechanical_neighs_idx
 * @return PetscErrorCode
 */
PetscErrorCode
DMSwarmRestoreParticleNeighbors(ParticleTopology *particle_topology,
                                IS mechanical_neighs_idx);

#endif /* NEIGHBORS_BLOBS_HPP */

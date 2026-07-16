/**
 * @file Blobs/Ghosts.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief Ghost blob creation / destruction for Blob3D.
 * @version 0.1
 * @date 2024-08-22
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef GHOST_BLOBS_HPP
#define GHOST_BLOBS_HPP

#include "Simulation.hpp"

/**
 * @brief Add ghost blobs to the simulation
 *
 * @param simulation
 * @param buffer_width
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmCreateGhostBlobs(Simulation &simulation,
                                       double buffer_width);

/**
 * @brief Remove ghost points
 *
 * @param simulation
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmDestroyGhostBlobs(Simulation &simulation);

#endif /* GHOST_BLOBS_HPP */

/**
 * @file Blobs/Ghosts.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2024-08-22
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "Macros.hpp"

#ifndef GHOST_BLOBS_HPP
#define GHOST_BLOBS_HPP

/**
 * @brief Add ghost atoms to the simulation
 *
 * @param Simulation
 * @param buffer_width
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmCreateGhostBlobs(DMD* Simulation, double buffer_width);

/**
 * @brief Remove ghost points
 *
 * @param Simulation
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmDestroyGhostBlobs(DMD* Simulation);

#endif /* GHOST_BLOBS_HPP */
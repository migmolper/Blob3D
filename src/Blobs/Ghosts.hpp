/**
 * @file Blobs/Ghosts.hpp
 */

#ifndef GHOST_BLOBS_HPP
#define GHOST_BLOBS_HPP

#include "Simulation.hpp"

PetscErrorCode DMSwarmCreateGhostBlobs(Simulation& simulation,
                                       double buffer_width);

PetscErrorCode DMSwarmDestroyGhostBlobs(Simulation& simulation);

#endif /* GHOST_BLOBS_HPP */

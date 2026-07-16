/**
 * @file Blobs/Topology.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2024-08-02
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "Macros.hpp"
#include "petscis.h"
#include <petscerror.h>

#ifndef TOPOLOGY_BLOBS_HPP
#define TOPOLOGY_BLOBS_HPP

/*******************************************************/

/**
 * @brief
 *
 * @param Simulation
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmGenerateBlobsTopology(DMD* Simulation, double r_cutoff);

/**
 * @brief
 *
 * @param Simulation
 * @param VOL_EXPANSION
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmRegenerateBlobsTopology(DMD* Simulation,
                                              double buffer_width);

/**
 * @brief
 *
 * @param Simulation
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmDestroyBlobsTopology(DMD* Simulation);

/*******************************************************/

#endif /* TOPOLOGY_BLOBS_HPP */
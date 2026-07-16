/**
 * @file Blobs/IO/vtk-output.hpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief
 * @version 0.1
 * @date 2022-06-22
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifdef USE_VTK

#ifndef OUTPUT_BLOBS_VTK_HPP
#define OUTPUT_BLOBS_VTK_HPP

#include "Macros.hpp"

/**
 * @brief
 *
 * @param step
 * @param adp
 * @return int
 */
PetscErrorCode DMSwarmBlobsViewVtk(DMD* Simulation, const std::string& filename);

#endif /* OUTPUT_BLOBS_VTK_HPP */

#endif /* OUTPUT_BLOBS_VTK_HPP */
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

#ifndef OUTPUT_BLOBS_VTK_HPP
#define OUTPUT_BLOBS_VTK_HPP

#include "Simulation.hpp"

/**
 * @brief
 *
 * @param step
 * @param adp
 * @return int
 */
PetscErrorCode DMSwarmBlobsViewVtk(Simulation& simulation,
                                   const std::string& filename);

#endif /* OUTPUT_BLOBS_VTK_HPP */

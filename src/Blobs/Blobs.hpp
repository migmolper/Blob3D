/**
 * @file Blobs/Blob.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BLOB_HPP
#define BLOB_HPP

#include "Macros.hpp"
#include <vector>

using namespace std;

/**
 * @brief
 *
 * @param Simulation
 * @param Simulation_file
 * @param IsingModel
 * @param r_cutoff_V
 * @return PetscErrorCode
 */
PetscErrorCode init_Blob_simulation(DMD* Simulation,               //!
                                    dump_file Simulation_file,     //!
                                    BackgroundMeshType mesh_type,  //!
                                    double r_cutoff_V);

/**
 * @brief
 *
 * @param Simulation
 */
PetscErrorCode destroy_Blob_simulation(DMD* Simulation);

#endif /* BLOB_HPP */
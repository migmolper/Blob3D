/**
 * @file Mesh/Create-Mesh.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef Create_Mesh_HPP
#define Create_Mesh_HPP

#include "Macros.hpp"
#include <vector>

using namespace std;

/**
 * @brief Create a Mesh object
 *
 * @param FE_Mesh
 * @param Simulation_file
 * @param r_cutoff_V Potential cutoff radius; stored on the mesh as buffer width
 * @return PetscErrorCode
 */
PetscErrorCode CreateMesh(DM* FE_Mesh, dump_file Simulation_file,
                          BackgroundMeshType type, double r_cutoff_V);

#endif /* Create_Mesh_HPP */
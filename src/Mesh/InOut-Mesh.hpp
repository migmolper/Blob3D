/**
 * @file Mesh/InOut-Mesh.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief Point-in-brick tests for mesh-aligned regions.
 * @version 0.1
 * @date 2025-01-16
 *
 * @copyright Copyright (c) 2025
 *
 */

#ifndef InOut_Mesh_HPP
#define InOut_Mesh_HPP

#include <Eigen/Dense>
#include <petscsnes.h>

/**
 * @brief Half-open brick test [lw, up) on each axis.
 *
 * Matches PETSc DMDA local partitioning and DMGetPhysicalLocalBoundingBox.
 * Use for periodic wrapping, ghost-buffer detection, and partition-aligned
 * regions.
 *
 * @param mean_q Point coordinates.
 * @param el_coords {x_lw, y_lw, z_lw, x_up, y_up, z_up}.
 */
bool In_Out_Mesh(const Eigen::Vector3d mean_q, const PetscScalar el_coords[]);

/**
 * @brief Closed brick test [lw, up] on each axis.
 *
 * Use for user-defined BC/fix regions where boundary atoms must be included.
 *
 * @param mean_q Point coordinates.
 * @param el_coords {x_lw, y_lw, z_lw, x_up, y_up, z_up}.
 */
bool In_Out_Mesh_Closed(const Eigen::Vector3d mean_q,
                        const PetscScalar el_coords[]);

#endif /* InOut_Mesh_HPP */

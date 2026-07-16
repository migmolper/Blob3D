/**
 * @file Mesh/InOut-Mesh.cpp
 * @author Miguel Molinos (@migmolper)
 * @brief Point-in-brick tests for mesh-aligned regions.
 * @version 0.1
 * @date 2025-01-16
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "Macros.hpp"
#include "petscdmplex.h"
#include <Eigen/Dense>
#include <petscdm.h>
#include <petscdmda.h>

/********************************************************************************/

/**
 * @brief Test membership in an axis-aligned brick.
 *
 * @param mean_q Point coordinates.
 * @param el_coords Brick bounds {x_lw, y_lw, z_lw, x_up, y_up, z_up}.
 * @param upper_inclusive If false, use half-open [lw, up) on each axis (PETSc
 *        DMDA partition convention). If true, use closed [lw, up].
 */
static bool In_Out_Brick(const Eigen::Vector3d mean_q,
                         const PetscScalar el_coords[], bool upper_inclusive) {

  const PetscReal xI_lw = el_coords[0];
  const PetscReal yI_lw = el_coords[1];
  const PetscReal zI_lw = el_coords[2];
  const PetscReal xI_up = el_coords[3];
  const PetscReal yI_up = el_coords[4];
  const PetscReal zI_up = el_coords[5];

  if (upper_inclusive) {
    return (mean_q(0) >= xI_lw) && (mean_q(0) <= xI_up) &&
           (mean_q(1) >= yI_lw) && (mean_q(1) <= yI_up) &&
           (mean_q(2) >= zI_lw) && (mean_q(2) <= zI_up);
  }

  return (mean_q(0) >= xI_lw) && (mean_q(0) < xI_up) &&
         (mean_q(1) >= yI_lw) && (mean_q(1) < yI_up) &&
         (mean_q(2) >= zI_lw) && (mean_q(2) < zI_up);
}

/********************************************************************************/

bool In_Out_Mesh(const Eigen::Vector3d mean_q, const PetscScalar el_coords[]) {
  return In_Out_Brick(mean_q, el_coords, false);
}

/********************************************************************************/

bool In_Out_Mesh_Closed(const Eigen::Vector3d mean_q,
                        const PetscScalar el_coords[]) {
  return In_Out_Brick(mean_q, el_coords, true);
}

/********************************************************************************/

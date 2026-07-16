/**
 * @file boundary_conditions.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "Macros.hpp"
#include "Mesh/Coordinates.hpp"
#include "Simulation.hpp"
#include <Eigen/Dense>

inline void MillerIndexToTuple(Miller_Index idx, int& dx, int& dy, int& dz) {
  switch (idx) {
    case idx_0_0_0:
      dx = 0;
      dy = 0;
      dz = 0;
      break;
    case idx_m1_m1_m1:
      dx = -1;
      dy = -1;
      dz = -1;
      break;
    case idx_m1_m1_0:
      dx = -1;
      dy = -1;
      dz = 0;
      break;
    case idx_m1_m1_p1:
      dx = -1;
      dy = -1;
      dz = +1;
      break;
    case idx_m1_0_m1:
      dx = -1;
      dy = 0;
      dz = -1;
      break;
    case idx_m1_0_0:
      dx = -1;
      dy = 0;
      dz = 0;
      break;
    case idx_m1_0_p1:
      dx = -1;
      dy = 0;
      dz = +1;
      break;
    case idx_m1_p1_m1:
      dx = -1;
      dy = +1;
      dz = -1;
      break;
    case idx_m1_p1_0:
      dx = -1;
      dy = +1;
      dz = 0;
      break;
    case idx_m1_p1_p1:
      dx = -1;
      dy = +1;
      dz = +1;
      break;
    case idx_0_m1_m1:
      dx = 0;
      dy = -1;
      dz = -1;
      break;
    case idx_0_m1_0:
      dx = 0;
      dy = -1;
      dz = 0;
      break;
    case idx_0_m1_p1:
      dx = 0;
      dy = -1;
      dz = +1;
      break;
    case idx_0_0_m1:
      dx = 0;
      dy = 0;
      dz = -1;
      break;
    case idx_0_0_p1:
      dx = 0;
      dy = 0;
      dz = +1;
      break;
    case idx_0_p1_m1:
      dx = 0;
      dy = +1;
      dz = -1;
      break;
    case idx_0_p1_0:
      dx = 0;
      dy = +1;
      dz = 0;
      break;
    case idx_0_p1_p1:
      dx = 0;
      dy = +1;
      dz = +1;
      break;
    case idx_p1_m1_m1:
      dx = +1;
      dy = -1;
      dz = -1;
      break;
    case idx_p1_m1_0:
      dx = +1;
      dy = -1;
      dz = 0;
      break;
    case idx_p1_m1_p1:
      dx = +1;
      dy = -1;
      dz = +1;
      break;
    case idx_p1_0_m1:
      dx = +1;
      dy = 0;
      dz = -1;
      break;
    case idx_p1_0_0:
      dx = +1;
      dy = 0;
      dz = 0;
      break;
    case idx_p1_0_p1:
      dx = +1;
      dy = 0;
      dz = +1;
      break;
    case idx_p1_p1_m1:
      dx = +1;
      dy = +1;
      dz = -1;
      break;
    case idx_p1_p1_0:
      dx = +1;
      dy = +1;
      dz = 0;
      break;
    case idx_p1_p1_p1:
      dx = +1;
      dy = +1;
      dz = +1;
      break;
  }
}

/*******************************************************/

inline bool IsDirectionActive(int dx, int dy, int dz, DMBoundaryType bx,
                              DMBoundaryType by, DMBoundaryType bz) {
  // dx direction
  if (dx != 0 && bx != DM_BOUNDARY_PERIODIC) return false;

  // dy direction
  if (dy != 0 && by != DM_BOUNDARY_PERIODIC) return false;

  // dz direction
  if (dz != 0 && bz != DM_BOUNDARY_PERIODIC) return false;

  return true;
}

/*******************************************************/

/**
 * @brief Get the mesh boundary condition object
 *
 * @param bcc Boundary condition in the (x, y, z) direction
 * @param da Mesh object
 * @return PetscErrorCode
 */
PetscErrorCode get_mesh_boundary_condition(DMBoundaryType* bcc, DM* da);

/**
 * @brief Enforce periodic boundary conditions over the ghost atoms
 *
 * @param Simulation Simulation object
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmEnforceBlobsPeriodic(Simulation& simulation,
                                           double r_cutoff);

/**
 * @brief Enforce periodic boundary conditions over the ghost atoms
 *
 * @param simulation
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmEnforceGhostBlobsPeriodic(Simulation& simulation);

/**
 * @brief
 *
 * @param mean_q
 * @param box_idx_ptr
 * @param background_mesh
 * @return PetscErrorCode
 */
PetscErrorCode VecEnforceGhostBlobsPeriodic(Vec mean_q,
                                            const PetscInt* box_idx_ptr,
                                            DM background_mesh);

PetscErrorCode DMSwarmFixMeanPositionBox(Simulation& simulation,
                                         PetscInt FixLabel,
                                         const PetscScalar box_coords[6]);

PetscErrorCode DMSwarmApplyDisplacement(Simulation& simulation,
                                        PetscInt FixLabel,
                                        PetscScalar displacement_x,
                                        PetscScalar displacement_y,
                                        PetscScalar displacement_z);

PetscErrorCode VecFixMeanPositionRHS(Vec RHS, Vec mean_q, Vec mean_q_ref,
                                     const PetscInt* idx_bcc_mean_q);

PetscErrorCode DMSwarmFixChemicalMultiplierBox(Simulation& simulation,
                                               const PetscScalar box_coords[6]);

PetscErrorCode DMSwarmFixThermalMultiplierBox(Simulation& simulation,
                                              const PetscScalar box_coords[6]);

#endif /* BOUNDARY_CONDITIONS_HPP */

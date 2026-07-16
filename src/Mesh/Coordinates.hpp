/**
 * @file Mesh/Coordinates.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief Reference ([-1,1]^3) vs physical coordinate utilities for DMDA/DMSwarm.
 * @version 0.1
 * @date 2026-06-22
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef COORDINATES_HPP
#define COORDINATES_HPP

#include "Macros.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <petscdm.h>

/**
 * Dual-coordinate layout used across SOLERA:
 *
 * - DMDA reference coordinates (`DMGetCoordinates`) stay fixed in [-1,1]^3 and
 *   drive PETSc partitioning, periodicity and `DMSwarmMigrate`.
 * - The global physical bounding box is stored on the DMDA (`_DMD_physical_bbox`)
 *   and maps reference space to physical units for BCs and PIC normalization.
 * - The neighbor-list buffer width is stored on the same DMDA container and set
 *   in CreateMesh from the potential cutoff radius.
 * - Swarm `mean-q` stores physical atom positions; `DMSwarmPICField_coor`
 *   stores the normalized positions required by the swarm PIC layer.
 */

/**
 * @brief Map a physical point to normalized [-1,1]^3 coordinates.
 *
 * @param gmin Global physical lower bound per direction
 * @param gmax Global physical upper bound per direction
 * @param real_x Physical x coordinate
 * @param real_y Physical y coordinate
 * @param real_z Physical z coordinate
 * @param norm_x Normalized x coordinate (output)
 * @param norm_y Normalized y coordinate (output)
 * @param norm_z Normalized z coordinate (output)
 */
void RealToNormalized(const PetscReal gmin[3], const PetscReal gmax[3],
                    PetscScalar real_x, PetscScalar real_y, PetscScalar real_z,
                    PetscScalar* norm_x, PetscScalar* norm_y,
                    PetscScalar* norm_z);

/**
 * @brief Map a normalized [-1,1]^3 point to physical coordinates.
 *
 * @param gmin Global physical lower bound per direction
 * @param gmax Global physical upper bound per direction
 * @param norm_x Normalized x coordinate
 * @param norm_y Normalized y coordinate
 * @param norm_z Normalized z coordinate
 * @param real_x Physical x coordinate (output)
 * @param real_y Physical y coordinate (output)
 * @param real_z Physical z coordinate (output)
 */
void NormalizedToReal(const PetscReal gmin[3], const PetscReal gmax[3],
                    PetscScalar norm_x, PetscScalar norm_y, PetscScalar norm_z,
                    PetscScalar* real_x, PetscScalar* real_y,
                    PetscScalar* real_z);

/**
 * @brief Global physical bounding box stored on the DMDA.
 *
 * @param dm DMDA background mesh
 * @param gmin Global physical lower bound (output)
 * @param gmax Global physical upper bound (output)
 * @param buffer_width Neighbor-list halo width in physical units (output, or NULL)
 * @return PetscErrorCode
 */
PetscErrorCode DMGetPhysicalBoundingBox(DM dm, PetscReal gmin[3],
                                        PetscReal gmax[3],
                                        PetscReal* buffer_width);

/**
 * @brief Local physical bounding box of the MPI partition.
 *
 * Affine map from the reference local bbox in [-1,1]^3 to physical space.
 *
 * @param dm DMDA background mesh
 * @param lmin Local physical lower bound (output)
 * @param lmax Local physical upper bound (output)
 * @return PetscErrorCode
 */
PetscErrorCode DMGetPhysicalLocalBoundingBox(DM dm, PetscReal lmin[3],
                                             PetscReal lmax[3]);

/**
 * @brief Initialize reference DMDA coordinates and store the physical bbox.
 *
 * Sets reference coordinates to [-1,1]^3 and records the input-file physical
 * bounding box on the DM for DMGetPhysicalBoundingBox.
 *
 * @param dm DMDA background mesh
 * @param gmin Global physical lower bound
 * @param gmax Global physical upper bound
 * @param buffer_width Neighbor-list halo width in physical units
 * @return PetscErrorCode
 */
PetscErrorCode DMCreatePhysicalCoordinates(DM dm, const PetscReal gmin[3],
                                           const PetscReal gmax[3],
                                           PetscReal buffer_width);

/**
 * @brief Store the neighbor-list buffer width on the DMDA.
 *
 * Typically set once in CreateMesh from the potential cutoff radius.
 *
 * @param dm DMDA background mesh
 * @param buffer_width Ghost/neighbor halo width in physical units
 * @return PetscErrorCode
 */
PetscErrorCode DMSetBufferWidth(DM dm, PetscReal buffer_width);

/**
 * @brief Read the neighbor-list buffer width stored on the DMDA.
 *
 * @param dm DMDA background mesh
 * @param buffer_width Buffer width in physical units (output)
 * @return PetscErrorCode
 */
PetscErrorCode DMGetBufferWidth(DM dm, PetscReal* buffer_width);

/**
 * @brief Scale the stored physical bbox on periodic directions.
 *
 * @param dm DMDA background mesh
 * @param F Deformation gradient (diagonal components used)
 * @return PetscErrorCode
 */
PetscErrorCode DMScalePhysicalCoordinates(DM dm, const Eigen::Matrix3d& F);

/**
 * @brief Build a diagonal deformation gradient from solver unknowns.
 *
 * Non-periodic directions are left at unity.
 */
inline Eigen::Matrix3d DMBuildDeformationGradient(const PetscScalar F_diag[3],
                                                   const DMBoundaryType bcc[3]) {
  Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
  for (unsigned int d = 0; d < NumberDimensions; d++) {
    if (bcc[d] == DM_BOUNDARY_PERIODIC) {
      F(d, d) = F_diag[d];
    }
  }
  return F;
}

/**
 * @brief Scale the stored physical bbox on the background mesh.
 *
 * After this call, `DMGetPhysicalBoundingBox` returns the deformed domain
 * limits (periodic directions scaled by the diagonal of F).
 *
 * @param Simulation Simulation object
 * @param F Deformation gradient (diagonal components used)
 * @return PetscErrorCode
 */
PetscErrorCode DMDApplyVolumetricExpansion(DMD* Simulation,
                                         const Eigen::Matrix3d& F);

/**
 * @brief Sync swarm PIC coordinates from physical `mean-q`.
 *
 * @param Simulation Simulation object
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmSyncCoorFromMeanQ(DMD* Simulation);

/**
 * @brief Sync swarm PIC coordinates from physical `mean-q`.
 *
 * @param atomistic_data DMSwarm object
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmSyncCoorFromMeanQ_DM(DM atomistic_data);

/**
 * @brief Scale physical `mean-q` with a diagonal deformation gradient.
 *
 * Only periodic directions are scaled. Used to commit the converged bulk
 * deformation gradient to the DMSwarm `mean-q` field after a barostat solve;
 * during the solve, trial positions live in a duplicated `Vec mean_q`.
 *
 * @param Simulation Simulation object
 * @param F Deformation gradient (diagonal components used)
 * @return PetscErrorCode
 */
PetscErrorCode DMSwarmScaleMeanQ(DMD* Simulation, const Eigen::Matrix3d& F);

#endif /* COORDINATES_HPP */

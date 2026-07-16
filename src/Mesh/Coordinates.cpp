/**
 * @file Mesh/Coordinates.cpp
 * @author Miguel Molinos (@migmolper)
 * @brief Reference ([-1,1]^3) vs physical coordinate utilities for
 * DMDA/DMSwarm.
 * @version 0.1
 * @date 2026-06-22
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "Mesh/Coordinates.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include <petscdmda.h>
#include <petscdmswarm.h>
#include <petscsys.h>

/************************************************************************/

void RealToNormalized(const PetscReal gmin[3], const PetscReal gmax[3],
                    PetscScalar real_x, PetscScalar real_y, PetscScalar real_z,
                    PetscScalar* norm_x, PetscScalar* norm_y,
                    PetscScalar* norm_z) {
  for (unsigned int d = 0; d < NumberDimensions; d++) {
    const PetscReal real[3] = {real_x, real_y, real_z};
    const PetscReal span = gmax[d] - gmin[d];
    PetscScalar normalized = 0.0;
    if (span > 0.0)
      normalized = 2.0 * (real[d] - gmin[d]) / span - 1.0;
    if (d == 0)
      *norm_x = normalized;
    else if (d == 1)
      *norm_y = normalized;
    else
      *norm_z = normalized;
  }
}

/************************************************************************/

void NormalizedToReal(const PetscReal gmin[3], const PetscReal gmax[3],
                    PetscScalar norm_x, PetscScalar norm_y, PetscScalar norm_z,
                    PetscScalar* real_x, PetscScalar* real_y,
                    PetscScalar* real_z) {
  const PetscScalar normalized[3] = {norm_x, norm_y, norm_z};
  PetscScalar* real[3] = {real_x, real_y, real_z};
  for (unsigned int d = 0; d < NumberDimensions; d++) {
    const PetscReal span = gmax[d] - gmin[d];
    *real[d] = gmin[d] + 0.5 * (normalized[d] + 1.0) * span;
  }
}

/************************************************************************/

/**
 * @brief Physical metadata stored on the DMDA via PetscObjectCompose.
 *
 * Attached to the DM as "_DMD_physical_bbox". Reference coordinates remain in
 * [-1,1]^3; this struct maps them to the physical brick used by BCs and PIC
 * normalization.
 */
typedef struct {
  PetscReal gmin[3];
  PetscReal gmax[3];
  PetscReal buffer_width;
  PetscBool has_buffer_width;
} DMPhysicalBBox;

#if PETSC_VERSION_LT(3, 25, 0)
static PetscErrorCode DMPhysicalBBoxDestroy(void** ctx) {
  PetscFunctionBeginUser;
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#else
static PetscErrorCode DMPhysicalBBoxDestroy(PetscCtxRt ctx) {
  PetscFunctionBeginUser;
  PetscCall(PetscFree(ctx));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

/************************************************************************/

static PetscErrorCode DMSetStoredPhysicalBoundingBox(DM dm,
                                                     const PetscReal gmin[3],
                                                     const PetscReal gmax[3],
                                                     PetscReal buffer_width) {
  DMPhysicalBBox* bbox = NULL;
  PetscContainer container = NULL;

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get or create the PetscContainer on the DM
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(PetscObjectQuery((PetscObject)dm, "_DMD_physical_bbox",
                             (PetscObject*)&container));
  if (container) {
    PetscCall(PetscContainerGetPointer(container, (void**)&bbox));
  } else {
    PetscCall(PetscNew(&bbox));
    bbox->has_buffer_width = PETSC_FALSE;
    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &container));
    PetscCall(PetscContainerSetPointer(container, bbox));
    PetscCall(PetscContainerSetCtxDestroy(container, DMPhysicalBBoxDestroy));
    PetscCall(PetscObjectCompose((PetscObject)dm, "_DMD_physical_bbox",
                                 (PetscObject)container));
    PetscCall(PetscObjectDereference((PetscObject)container));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Store physical bbox and neighbor-list buffer width
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (unsigned int d = 0; d < NumberDimensions; d++) {
    bbox->gmin[d] = gmin[d];
    bbox->gmax[d] = gmax[d];
  }
  bbox->buffer_width = buffer_width;
  bbox->has_buffer_width = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

static PetscErrorCode DMGetStoredPhysicalBoundingBox(DM dm, PetscReal gmin[3],
                                                     PetscReal gmax[3],
                                                     PetscReal* buffer_width,
                                                     PetscBool* has_bbox) {
  PetscContainer container = NULL;
  DMPhysicalBBox* bbox = NULL;

  PetscFunctionBeginUser;
  *has_bbox = PETSC_FALSE;
  PetscCall(PetscObjectQuery((PetscObject)dm, "_DMD_physical_bbox",
                             (PetscObject*)&container));
  if (!container) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscContainerGetPointer(container, (void**)&bbox));
  for (unsigned int d = 0; d < NumberDimensions; d++) {
    gmin[d] = bbox->gmin[d];
    gmax[d] = bbox->gmax[d];
  }
  if (buffer_width) {
    PetscCheck(bbox->has_buffer_width, PetscObjectComm((PetscObject)dm),
               PETSC_ERR_ARG_WRONGSTATE,
               "Buffer width not set; call DMCreatePhysicalCoordinates first");
    *buffer_width = bbox->buffer_width;
  }
  *has_bbox = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSetBufferWidth(DM dm, PetscReal buffer_width) {
  PetscReal gmin[3], gmax[3], current_buffer_width;

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Read current bbox and overwrite the stored buffer width
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(
      DMGetPhysicalBoundingBox(dm, gmin, gmax, &current_buffer_width));
  PetscCall(DMSetStoredPhysicalBoundingBox(dm, gmin, gmax, buffer_width));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMGetBufferWidth(DM dm, PetscReal* buffer_width) {
  PetscReal gmin[3], gmax[3];

  PetscFunctionBeginUser;
  PetscCheck(buffer_width, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_NULL,
             "buffer_width must not be NULL");
  PetscCall(DMGetPhysicalBoundingBox(dm, gmin, gmax, buffer_width));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMGetPhysicalBoundingBox(DM dm, PetscReal gmin[3],
                                        PetscReal gmax[3],
                                        PetscReal* buffer_width) {
  PetscBool has_bbox = PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCall(DMGetStoredPhysicalBoundingBox(dm, gmin, gmax, buffer_width,
                                           &has_bbox));
  PetscCheck(has_bbox, PetscObjectComm((PetscObject)dm),
             PETSC_ERR_ARG_WRONGSTATE,
             "Physical bounding box not set; call DMCreatePhysicalCoordinates "
             "first");
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMGetPhysicalLocalBoundingBox(DM dm, PetscReal lmin[3],
                                             PetscReal lmax[3]) {
  PetscReal lmin_ref[3], lmax_ref[3];
  PetscReal gmin[3], gmax[3];

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Map the reference local bbox in [-1,1]^3 to physical space
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMGetLocalBoundingBox(dm, lmin_ref, lmax_ref));
  PetscCall(DMGetPhysicalBoundingBox(dm, gmin, gmax, NULL));
  for (unsigned int d = 0; d < NumberDimensions; d++) {
    const PetscReal span = gmax[d] - gmin[d];
    lmin[d] = gmin[d] + 0.5 * (lmin_ref[d] + 1.0) * span;
    lmax[d] = gmin[d] + 0.5 * (lmax_ref[d] + 1.0) * span;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMCreatePhysicalCoordinates(DM dm, const PetscReal gmin[3],
                                           const PetscReal gmax[3],
                                           PetscReal buffer_width) {
  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set reference DMDA coordinates and store the physical bbox
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMDASetUniformCoordinates(dm, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0));
  PetscCall(DMSetStoredPhysicalBoundingBox(dm, gmin, gmax, buffer_width));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMScalePhysicalCoordinates(DM dm, const Eigen::Matrix3d& F) {
  PetscReal gmin[3], gmax[3], buffer_width;
  PetscReal new_gmin[3], new_gmax[3];
  DMBoundaryType bx, by, bz;

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the global physical bbox and mesh boundary conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMGetPhysicalBoundingBox(dm, gmin, gmax, &buffer_width));
  PetscCall(DMDAGetInfo(dm, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                        NULL, &bx, &by, &bz, NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Scale periodic directions only
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (unsigned int d = 0; d < NumberDimensions; d++) {
    new_gmin[d] = gmin[d];
    new_gmax[d] = gmax[d];
  }
  if (bx == DM_BOUNDARY_PERIODIC) {
    new_gmin[0] *= F(0, 0);
    new_gmax[0] *= F(0, 0);
  }
  if (by == DM_BOUNDARY_PERIODIC) {
    new_gmin[1] *= F(1, 1);
    new_gmax[1] *= F(1, 1);
  }
  if (bz == DM_BOUNDARY_PERIODIC) {
    new_gmin[2] *= F(2, 2);
    new_gmax[2] *= F(2, 2);
  }

  PetscCall(DMSetStoredPhysicalBoundingBox(dm, new_gmin, new_gmax,
                                           buffer_width));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMDApplyVolumetricExpansion(DMD* Simulation,
                                           const Eigen::Matrix3d& F) {
  DM background_mesh;

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Scale the stored physical bbox on the background mesh
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmGetCellDM(Simulation->atomistic_data, &background_mesh));
  PetscCall(DMScalePhysicalCoordinates(background_mesh, F));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmSyncCoorFromMeanQ_DM(DM atomistic_data) {
  PetscFunctionBeginUser;

  DM background_mesh;
  PetscReal gmin[3], gmax[3];
  PetscScalar* mean_q_ptr;
  PetscScalar* PIC_coor_ptr;
  PetscInt n_local;
  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Physical bbox
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmGetCellDM(atomistic_data, &background_mesh));
  PetscCall(DMGetPhysicalBoundingBox(background_mesh, gmin, gmax, NULL));
  PetscCall(DMSwarmGetLocalSize(atomistic_data, &n_local));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get swarm fields
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmGetField(atomistic_data, "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  PetscCall(DMSwarmGetField(atomistic_data, DMSwarmPICField_coor, NULL, NULL,
                            (void**)&PIC_coor_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    mean-q (physical) -> DMSwarmPICField_coor (normalized)
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (PetscInt site_i = 0; site_i < n_local; site_i++) {
    RealToNormalized(gmin, gmax,                       //!
                     mean_q_ptr[site_i * dim + 0],     //!
                     mean_q_ptr[site_i * dim + 1],     //!
                     mean_q_ptr[site_i * dim + 2],     //!
                     &PIC_coor_ptr[site_i * dim + 0],  //!
                     &PIC_coor_ptr[site_i * dim + 1],  //!
                     &PIC_coor_ptr[site_i * dim + 2]);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore swarm fields
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(atomistic_data, "mean-q", NULL, NULL,
                                (void**)&mean_q_ptr));
  PetscCall(DMSwarmRestoreField(atomistic_data, DMSwarmPICField_coor, NULL,
                                NULL, (void**)&PIC_coor_ptr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmSyncCoorFromMeanQ(DMD* Simulation) {
  PetscFunctionBeginUser;
  PetscCall(DMSwarmSyncCoorFromMeanQ_DM(Simulation->atomistic_data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode DMSwarmScaleMeanQ(DMD* Simulation, const Eigen::Matrix3d& F) {
  DM background_mesh;
  DMBoundaryType bcc[3] = {DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DM_BOUNDARY_NONE};
  PetscScalar* mean_q_ptr;
  PetscInt n_local;
  unsigned int dim = NumberDimensions;

  PetscFunctionBeginUser;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get mesh boundary conditions and mean-q field
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmGetCellDM(Simulation->atomistic_data, &background_mesh));
  PetscCall(get_mesh_boundary_condition(bcc, &background_mesh));
  const DMBoundaryType bx = bcc[0];
  const DMBoundaryType by = bcc[1];
  const DMBoundaryType bz = bcc[2];
  PetscCall(DMSwarmGetLocalSize(Simulation->atomistic_data, &n_local));
  PetscCall(DMSwarmGetField(Simulation->atomistic_data, "mean-q", NULL, NULL,
                            (void**)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Scale periodic directions only
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (PetscInt site_i = 0; site_i < n_local; site_i++) {
    if (bx == DM_BOUNDARY_PERIODIC) mean_q_ptr[site_i * dim + 0] *= F(0, 0);
    if (by == DM_BOUNDARY_PERIODIC) mean_q_ptr[site_i * dim + 1] *= F(1, 1);
    if (bz == DM_BOUNDARY_PERIODIC) mean_q_ptr[site_i * dim + 2] *= F(2, 2);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore mean-q data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(DMSwarmRestoreField(Simulation->atomistic_data, "mean-q", NULL,
                                NULL, (void**)&mean_q_ptr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

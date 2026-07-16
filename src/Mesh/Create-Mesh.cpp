/**
 * @file Mesh/Create-Mesh.cpp
 * @author Miguel Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

 #include <petscis.h>
 #include <petscsystypes.h>
 #if __APPLE__
 #include <malloc/_malloc.h>
 #endif
 #ifdef USE_MPI
 #include <mpi.h>
 #endif
 #include "Mesh/Coordinates.hpp"
 #include <Eigen/Dense>
 #include <iostream>
 #include <math.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 
 extern PetscMPIInt size_MPI;
 extern PetscMPIInt rank_MPI;
 
 extern PetscInt ndiv_mesh_X;
 extern PetscInt ndiv_mesh_Y;
 extern PetscInt ndiv_mesh_Z;
 
 extern PetscInt size_MPI_X;
 extern PetscInt size_MPI_Y;
 extern PetscInt size_MPI_Z;
 
 using namespace std;
 
 /*******************************************************/
 
 static PetscErrorCode CreateMesh_DMDA(DM* FE_Mesh,                //!
                                       dump_file Simulation_file,  //!
                                       double r_cutoff_V           //!
 ) {
 
   PetscReal gmin[3] = {0, 0., 0.}, gmax[3] = {0, 0., 0.};
 
   PetscFunctionBeginUser;
 
   unsigned int dim = NumberDimensions;
   PetscInt overlap = 1;
 
   double box_x_min = Simulation_file.box_x_min;
   double box_x_max = Simulation_file.box_x_max;
   double box_y_min = Simulation_file.box_y_min;
   double box_y_max = Simulation_file.box_y_max;
   double box_z_min = Simulation_file.box_z_min;
   double box_z_max = Simulation_file.box_z_max;
 
   PetscCall(DMDACreate(PETSC_COMM_WORLD, FE_Mesh));
   PetscCall(DMSetDimension(*FE_Mesh, dim));
   PetscCall(DMDASetSizes(*FE_Mesh, ndiv_mesh_X, ndiv_mesh_Y, ndiv_mesh_Z));
   PetscCall(DMDASetNumProcs(*FE_Mesh, size_MPI_X, size_MPI_Y, size_MPI_Z));
   PetscCall(DMDASetBoundaryType(*FE_Mesh, Simulation_file.bx,
                                 Simulation_file.by, Simulation_file.bz));
   PetscCall(DMDASetDof(*FE_Mesh, 1));
   PetscCall(DMDASetStencilType(*FE_Mesh, DMDA_STENCIL_BOX));
   PetscCall(DMDASetStencilWidth(*FE_Mesh, overlap));
   PetscCall(DMDASetElementType(*FE_Mesh, DMDA_ELEMENT_Q1));
   PetscCall(DMSetFromOptions(*FE_Mesh));
   PetscCall(DMSetUp(*FE_Mesh));
 
   //! Periodic directions need >= 2 MPI ranks: ghost atoms are exchanged across
   //! inter-rank faces; a single rank cannot provide a halo across the wrap.
   {
     PetscInt nprocs_x = 0, nprocs_y = 0, nprocs_z = 0;
     DMBoundaryType bx = Simulation_file.bx;
     DMBoundaryType by = Simulation_file.by;
     DMBoundaryType bz = Simulation_file.bz;
 
     PetscCall(DMDAGetInfo(*FE_Mesh, NULL, NULL, NULL, NULL, &nprocs_x,
                           &nprocs_y, &nprocs_z, NULL, NULL, &bx, &by, &bz,
                           NULL));
 
     if (bx == DM_BOUNDARY_PERIODIC) {
       PetscCheck(
           nprocs_x >= 2, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
           "Periodic X requires at least 2 MPI ranks in X (got %" PetscInt_FMT
           ")",
           nprocs_x);
     }
     if (by == DM_BOUNDARY_PERIODIC) {
       PetscCheck(
           nprocs_y >= 2, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
           "Periodic Y requires at least 2 MPI ranks in Y (got %" PetscInt_FMT
           ")",
           nprocs_y);
     }
     if (bz == DM_BOUNDARY_PERIODIC) {
       PetscCheck(
           nprocs_z >= 2, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
           "Periodic Z requires at least 2 MPI ranks in Z (got %" PetscInt_FMT
           ")",
           nprocs_z);
     }
   }
 
   {
     PetscInt M = 0, N = 0, P = 0;
     PetscInt nprocs_x = 0, nprocs_y = 0, nprocs_z = 0;
     PetscCall(DMDAGetInfo(*FE_Mesh, NULL, &M, &N, &P, &nprocs_x, &nprocs_y,
                           &nprocs_z, NULL, NULL, NULL, NULL, NULL, NULL));
     PetscCheck(M >= nprocs_x, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "DMDA X grid points (%" PetscInt_FMT
                ") must be >= MPI ranks in X (%" PetscInt_FMT ")",
                M, nprocs_x);
     PetscCheck(N >= nprocs_y, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "DMDA Y grid points (%" PetscInt_FMT
                ") must be >= MPI ranks in Y (%" PetscInt_FMT ")",
                N, nprocs_y);
     PetscCheck(P >= nprocs_z, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "DMDA Z grid points (%" PetscInt_FMT
                ") must be >= MPI ranks in Z (%" PetscInt_FMT ")",
                P, nprocs_z);
   }
 
   const PetscReal physical_gmin[3] = {box_x_min, box_y_min, box_z_min};
   const PetscReal physical_gmax[3] = {box_x_max, box_y_max, box_z_max};
   PetscCall(DMCreatePhysicalCoordinates(*FE_Mesh, physical_gmin, physical_gmax,
                                         r_cutoff_V));
 
   {
     DMBoundaryType bx = Simulation_file.bx;
     DMBoundaryType by = Simulation_file.by;
     DMBoundaryType bz = Simulation_file.bz;
 
     if (bx == DM_BOUNDARY_PERIODIC || by == DM_BOUNDARY_PERIODIC ||
         bz == DM_BOUNDARY_PERIODIC) {
       //! Reference DMDA coordinates live in [-1,1]^3; periodicity must match.
       PetscReal maxCell[3] = {0., 0., 0.};
       PetscReal Lstart[3] = {-1.0, -1.0, -1.0};
       PetscReal L[3] = {0., 0., 0.};
       if (bx == DM_BOUNDARY_PERIODIC) {
         maxCell[0] = 2.0 / ndiv_mesh_X;
         L[0] = 2.0;
       }
       if (by == DM_BOUNDARY_PERIODIC) {
         maxCell[1] = 2.0 / ndiv_mesh_Y;
         L[1] = 2.0;
       }
       if (bz == DM_BOUNDARY_PERIODIC) {
         maxCell[2] = 2.0 / ndiv_mesh_Z;
         L[2] = 2.0;
       }
       PetscCall(DMSetPeriodicity(*FE_Mesh, maxCell, Lstart, L));
     }
   }
 
   //! @brief Check if the local dimensions are smaller than twice the cutoff
   //! radius
   PetscReal ll_vertex[3], ur_vertex[3];
   PetscCall(DMGetPhysicalLocalBoundingBox(*FE_Mesh, ll_vertex, ur_vertex));
 
   //! @brief check if the max and min of the local bounding box are smaller than
   //! twice the cutoff radius
   PetscCall(DMGetPhysicalBoundingBox(*FE_Mesh, gmin, gmax, NULL));
 
   if (fabs(box_x_max - gmax[0]) > 1e-6) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "Does not match the box_x_max (%f) match the global bbox (%f)?",
                          box_x_max, gmax[0]));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
 
   if (fabs(box_y_max - gmax[1]) > 1e-6) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "Does not match the box_y_max (%f) match the global bbox (%f)?",
                          box_y_max, gmax[1]));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
   if (fabs(box_z_max - gmax[2]) > 1e-6) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "Does not match the box_z_max (%f) match the global bbox (%f)?",
                          box_z_max, gmax[2]));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
   if (fabs(box_x_min - gmin[0]) > 1e-6) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "Does not match the box_x_min (%f) match the global bbox (%f)?",
                          box_x_min, gmin[0]));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
   if (fabs(box_y_min - gmin[1]) > 1e-6) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "Does not match the box_y_min (%f) match the global bbox (%f)?",
                          box_y_min, gmin[1]));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
   if (fabs(box_z_min - gmin[2]) > 1e-6) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "Does not match the box_z_min (%f) match the global bbox (%f)?",
                          box_z_min, gmin[2]));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
 
   PetscReal local_dx = ur_vertex[0] - ll_vertex[0];
   if (local_dx < 2 * r_cutoff_V) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "The local dimension X (%f) is smaller than %f",
                          local_dx, 2 * r_cutoff_V));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
 
   PetscReal local_dy = ur_vertex[1] - ll_vertex[1];
   if (local_dy < 2 * r_cutoff_V) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "The local dimension Y (%f) is smaller than %f",
                          local_dy, 2 * r_cutoff_V));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
 
   PetscReal local_dz = ur_vertex[2] - ll_vertex[2];
   if (local_dz < 2 * r_cutoff_V) {
     PetscCall(PetscError(PETSC_COMM_WORLD, __LINE__, "init_DMD_simulation",
                          __FILE__, PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                          "The local dimension Z (%f) is smaller than %f",
                          local_dz, 2 * r_cutoff_V));
     PetscFunctionReturn(PETSC_ERR_RETURN);
   }
 
   PetscFunctionReturn(PETSC_SUCCESS);
 }
 
 /*******************************************************/
 
 static PetscErrorCode CreateMesh_DMPLEX(DM* FE_Mesh,                //!
                                         dump_file Simulation_file,  //!
                                         double r_cutoff_V           //!
 ) {
 
   PetscFunctionBeginUser;
 
   unsigned int dim = NumberDimensions;
   PetscInt overlap = 1;
   PetscInt dof = 1;
   PetscInt localizationHeight = 0;
   PetscBool sparseLocalize = PETSC_TRUE;
   PetscBool simplex = PETSC_TRUE;
   PetscBool interpolate = PETSC_FALSE;
   const PetscInt faces[3] = {ndiv_mesh_X, ndiv_mesh_Y, ndiv_mesh_Z};
   const DMBoundaryType periodicity[3] = {Simulation_file.bx, Simulation_file.by,
                                          Simulation_file.bz};
 
   double box_ll[3] = {Simulation_file.box_x_min, Simulation_file.box_y_min,
                       Simulation_file.box_z_min};
   double box_ur[3] = {Simulation_file.box_x_max, Simulation_file.box_y_max,
                       Simulation_file.box_z_max};
 
   PetscCall(DMCreate(PETSC_COMM_WORLD, FE_Mesh));
   PetscCall(DMSetType(*FE_Mesh, DMPLEX));
   PetscCall(DMSetFromOptions(*FE_Mesh));
   PetscCall(DMViewFromOptions(*FE_Mesh, NULL, "-dm_view"));
 
   //  PetscCall(DMPlexCreateBoxMesh(
   //      PETSC_COMM_WORLD, dim, simplex, faces, box_ll,
   //      box_ur, periodicity, interpolate, localizationHeight,
   //      sparseLocalize, FE_Mesh));
 
   PetscFunctionReturn(PETSC_SUCCESS);
 }
 
 /*******************************************************/
 
 /**
  * @brief Create a Mesh object
  *
  * @param FE_Mesh
  * @param Simulation_file
  * @return PetscErrorCode
  */
PetscErrorCode CreateMesh(DM* FE_Mesh, dump_file Simulation_file,
                          BackgroundMeshType type, double r_cutoff_V) {

  PetscFunctionBeginUser;

  switch (type) {
    case BackgroundMeshType::DMDA_mesh:
      PetscCall(CreateMesh_DMDA(FE_Mesh, Simulation_file, r_cutoff_V));
      break;
    case BackgroundMeshType::DMPLEX_mesh:
      /* DMPLEX box-mesh path is not fully wired yet; use DMDA. */
      PetscCall(CreateMesh_DMDA(FE_Mesh, Simulation_file, r_cutoff_V));
      break;
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}
 
 /*******************************************************/
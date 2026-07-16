
/**
 * @file Mesh/test/test-Boundary-Conditions.cpp
 * @author Miguel Molinos (migmolper)
 * @brief
 * @version 0.1
 * @date 2025-02-04
 *
 * @copyright Copyright (c) 2025
 *
 */

 #include "Blobs/Blobs.hpp"
 #include "Blobs/Neighbors.hpp"
 #include "Blobs/Topology.hpp"
 #include "Blobs/IO/dump-inputs.hpp"
 #include "Blobs/IO/vtk-outputs.hpp"
 #include "Kinetics/Mass-Transport-JKO-TAO.hpp"
 #include "Macros.hpp"
 #include "Mesh/Boundary-Conditions.hpp"
 #include "Potentials/Advection-Diff-OpenMP.hpp"
 #include "catch.hpp"
 #include <Eigen/Dense>
 #include <cmath>
 #include <iostream>  //std::cout//std::cin
 #include <limits>
 #include <math.h>
 #include <memory>
 #include <stdio.h>
 #include <stdlib.h>
 
 extern PetscMPIInt size_MPI;
 extern PetscMPIInt rank_MPI;
 
 extern PetscInt ndiv_mesh_X;
 extern PetscInt ndiv_mesh_Y;
 extern PetscInt ndiv_mesh_Z;
 
 extern PetscInt size_MPI_X;
 extern PetscInt size_MPI_Y;
 extern PetscInt size_MPI_Z;
 
 extern DiffusivePotential Potential_AD;
 
 /********************************************************************************/
 
 TEST_CASE("Test diffusion JKO algorithm. Specimen: Mg-hcp-sphere-R-20",
           "[Mass-Trasport-Diffusion][Mg-hcp-sphere-R-20]") {
 
   if (rank_MPI == 0)
     std::cout << "[Mass-Trasport-Diffusion][Mg-hcp-sphere-R-20]" << std::endl;
 
   unsigned int dim = NumberDimensions;
   const double Tolerance = 1.E-8;
   PetscScalar Delta_r = 6.0;
   PetscScalar kappa = 0.1;    //!
   PetscScalar Delta_t = 1;    //!
   PetscScalar m_p = 0.01;     //!
   PetscScalar rho_ref = 0.2;  //!
   PetscErrorCode ierr;
   const char OutputFolder[MAXC] = "Kinetics/outputs";
   char SimulationFile[MAXC] = "Kinetics/inputs/Mg-hcp-sphere-R-20.dump";
   char OutputFile[MAXC];
 
   ndiv_mesh_X = 20;
   ndiv_mesh_Y = 20;
   ndiv_mesh_Z = 20;
 
   size_MPI_X = 2;
   size_MPI_Y = 2;
   size_MPI_Z = 2;
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Command line options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_gatol", "1.e-6");
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_gttol", "1.e-9");
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_max_it", "100");
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_type", "cg");
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_cg_type", "prp");
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_max_funcs", "100");
   //  PetscOptionsSetValue(NULL, "-minJKO_dx_tao_monitor_globalization", "");
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_view", "");
   PetscOptionsSetValue(NULL, "-minJKO_dx_tao_converged_reason", "");
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Read information from dump file
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   dump_file Simulation_dump_data = DMSwarmBlobsReadDump(SimulationFile);
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize atomistic simulation
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   DMD Simulation;
   REQUIRE_NOTHROW(init_Blob_simulation(&Simulation, Simulation_dump_data,
                                        BackgroundMeshType::DMPLEX_mesh,
                                        r_cutoff_ADP_MgHx));
   REQUIRE_NOTHROW(DMSwarmSetMigrateType(Simulation.atomistic_data,
                                         DMSWARM_MIGRATE_DMCELLNSCATTER));
   REQUIRE_NOTHROW(DMSwarmMigrate(Simulation.atomistic_data, PETSC_TRUE));
 
   //! Set potential parameters
   Potential_AD.C = 0.0;
   Potential_AD.kappa = kappa;
   Potential_AD.rho_ref = rho_ref;
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free dump data
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   DMSwarmBlobsFreeDump(&Simulation_dump_data);
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Output data
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   REQUIRE_NOTHROW(DMView(Simulation.atomistic_data, PETSC_VIEWER_STDOUT_WORLD));
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set beta, mass and mean momentum value
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   PetscScalar beta_value = 2.0 / DSQR(Delta_r);
   Vec beta;
   REQUIRE_NOTHROW(DMSwarmCreateGlobalVectorFromField(Simulation.atomistic_data,
                                                      "beta", &beta));
   REQUIRE_NOTHROW(VecSet(beta, beta_value));
   REQUIRE_NOTHROW(DMSwarmDestroyGlobalVectorFromField(Simulation.atomistic_data,
                                                       "beta", &beta));
 
   Vec mass;
   REQUIRE_NOTHROW(DMSwarmCreateGlobalVectorFromField(Simulation.atomistic_data,
                                                      "mass", &mass));
   REQUIRE_NOTHROW(VecSet(mass, m_p));
   REQUIRE_NOTHROW(DMSwarmDestroyGlobalVectorFromField(Simulation.atomistic_data,
                                                       "mass", &mass));
 
   Vec P;
   REQUIRE_NOTHROW(DMSwarmCreateGlobalVectorFromField(Simulation.atomistic_data,
                                                      "mean-p", &P));
   REQUIRE_NOTHROW(VecSet(P, 0.0));
   REQUIRE_NOTHROW(DMSwarmDestroyGlobalVectorFromField(Simulation.atomistic_data,
                                                       "mean-p", &P));
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create ghost atoms and topology
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   REQUIRE_NOTHROW(DMSwarmGenerateBlobsTopology(&Simulation, Delta_r));
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize  equations
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   governing_equations system_equations = Advection_Diff_constructor();
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Print background mesh
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   DM background_mesh;
   REQUIRE_NOTHROW(
       DMSwarmGetCellDM(Simulation.atomistic_data, &background_mesh));
 
   PetscViewer viewer;
   REQUIRE_NOTHROW(PetscViewerCreate(PETSC_COMM_WORLD, &viewer));
   REQUIRE_NOTHROW(PetscViewerSetType(viewer, PETSCVIEWERVTK));
   REQUIRE_NOTHROW(PetscViewerFileSetMode(viewer, FILE_MODE_WRITE));
   snprintf(OutputFile, MAXC, "%s/Mg-hcp-sphere-R-20-mesh.vtu", OutputFolder);
   REQUIRE_NOTHROW(PetscViewerFileSetName(viewer, OutputFile));
   REQUIRE_NOTHROW(DMView(background_mesh, viewer));
   REQUIRE_NOTHROW(PetscViewerDestroy(&viewer));
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Outputs initial configuration
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   snprintf(OutputFile, MAXC, "%s/Mg-hcp-sphere-R-20-%i.vtp", OutputFolder, 0);
   REQUIRE_NOTHROW(DMSwarmBlobsViewVtk(&Simulation, OutputFile));
 
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Update mean position
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   for (unsigned int i = 1; i < 50; i++) {
 
     //! Run advection-diffusion blob solver
     REQUIRE_NOTHROW(Mass_Trasport_Advection_Diffusion(Delta_t, &Simulation,
                                                       system_equations));
 
     //! Update simulation topology
     REQUIRE_NOTHROW(DMSwarmRegenerateBlobsTopology(&Simulation, Delta_r,
                                                    PETSC_TRUE, PETSC_FALSE));
 
     //! Output
     snprintf(OutputFile, MAXC, "%s/Mg-hcp-sphere-R-20-%i.vtp", OutputFolder, i);
     REQUIRE_NOTHROW(DMSwarmBlobsViewVtk(&Simulation, OutputFile));
   }
 
   //! Destroy list of neighbors and other important information
   REQUIRE_NOTHROW(DMSwarmDestroyBlobsTopology(&Simulation));
 
   //! @brief Destroy atomistic
   REQUIRE_NOTHROW(destroy_Blob_simulation(&Simulation));
 }
 
 /********************************************************************************/
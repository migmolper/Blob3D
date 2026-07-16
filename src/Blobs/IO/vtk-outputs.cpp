/**
 * @file Blobs/IO/vtk-output.cpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief
 * @version 0.1
 * @date 2022-06-22
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Blobs/Blobs.hpp"
#include "Macros.hpp"
#include <Eigen/Dense>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>  //std::cout//std::cin
#include <math.h>
#include <sstream>  //std::istringstream
#include <string>   //std::string
#include <vector>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTensorGlyph.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>

using namespace std;
extern char OutputFolder[MAXC];

/********************************************************************************/

PetscErrorCode DMSwarmBlobsViewVtk(Simulation& simulation, const std::string& filename) {

  PetscFunctionBeginUser;
  MPI_Comm comm = PETSC_COMM_WORLD;

  PetscMPIInt rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  PetscInt dim = NumberDimensions;
  PetscInt totalParticles;
  PetscInt numParticlesLocal;

  //! Get the DMSwarm object
  DM atomistic_data = simulation.dm();

  // Get local number of particles and their coordinates
  PetscCall(DMSwarmGetSize(atomistic_data, &totalParticles));
  PetscCall(DMSwarmGetLocalSize(atomistic_data, &numParticlesLocal));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the global (x,y,z) indices of the lower left corner and size of the
    local region, excluding ghost points.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscReal lmin[3], lmax[3];
  DM background_mesh;
  PetscCall(DMSwarmGetCellDM(atomistic_data, &background_mesh));
  PetscCall(DMGetBoundingBox(background_mesh, lmin, lmax));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get index of the particles
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt* idx_ptr;
  PetscCall(
      DMSwarmGetField(atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Boolean variable for ghost particles
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt* ghost_ptr;
  PetscCall(
      DMSwarmGetField(atomistic_data, "ghost", NULL, NULL, (void**)&ghost_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the mean position
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const PetscReal* mean_q_ptr = nullptr;
  PetscCall(DMSwarmGetField(atomistic_data, DMSwarmPICField_coor, NULL, NULL,
                            (void**)&mean_q_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the molar fraction
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar* rho_ptr;
  PetscCall(DMSwarmGetField(atomistic_data, "rho", NULL, NULL,
                            (void**)&rho_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get the thermal Lagrange Multiplier (beta)
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar* beta_ptr;
  PetscCall(
      DMSwarmGetField(atomistic_data, "beta", NULL, NULL, (void**)&beta_ptr));

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Index for the thermal bcc
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt* idx_beta_bcc_ptr;
  PetscCall(DMSwarmGetField(atomistic_data, "idx-bcc-beta", NULL, NULL,
                            (void**)&idx_beta_bcc_ptr));

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Index for the MPI rank
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt* rank_ptr;
  PetscCall(DMSwarmGetField(atomistic_data, DMSwarmField_rank, NULL, NULL,
                            (void**)&rank_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Calculate the local number of ghost particles
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt numParticlesGhostLocal = 0;
  for (PetscInt i = 0; i < numParticlesLocal; ++i) {
    if (ghost_ptr[i] == 1) {
      numParticlesGhostLocal++;
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Build local buffers (without ghosts)
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt numSitesLocal = numParticlesLocal - numParticlesGhostLocal;

  std::vector<double> local_pos;            // dim * numSitesLocal
  std::vector<double> local_xi;             // 1 * numSitesLocal
  std::vector<double> local_beta;           // 1 * numSitesLocal
  std::vector<int> local_idx_beta_bcc;      // 1 * numSitesLocal
  std::vector<int> local_rank;              // 1 * numSitesLocal

  local_pos.reserve(dim * numSitesLocal);
  local_xi.reserve(numSitesLocal);
  local_beta.reserve(numSitesLocal);
  local_idx_beta_bcc.reserve(numSitesLocal);
  local_rank.reserve(numSitesLocal);
  
  for (PetscInt i = 0; i < numSitesLocal; ++i) {

    // positions
    local_pos.push_back(mean_q_ptr[i * dim + 0]);
    local_pos.push_back(mean_q_ptr[i * dim + 1]);
    local_pos.push_back(mean_q_ptr[i * dim + 2]);

    // scalar fields
    local_xi.push_back(PetscRealPart(rho_ptr[i]));
    local_beta.push_back(PetscRealPart(beta_ptr[i]));
    local_idx_beta_bcc.push_back(static_cast<int>(idx_beta_bcc_ptr[i]));
    local_rank.push_back(static_cast<int>(rank_ptr[i]));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Gather all sites on rank 0 (MPI_Gatherv)
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  std::vector<int> recvcounts_sites, displs_sites;
  int ns_local_int = static_cast<int>(numSitesLocal);

  if (rank == 0) {
    recvcounts_sites.resize(size);
  }

  MPI_Gather(&ns_local_int, 1, MPI_INT,
             rank == 0 ? recvcounts_sites.data() : nullptr, 1, MPI_INT, 0,
             comm);

  PetscInt numSitesGlobal = 0;
  if (rank == 0) {
    displs_sites.resize(size);
    displs_sites[0] = 0;
    numSitesGlobal = recvcounts_sites[0];
    for (int p = 1; p < size; ++p) {
      displs_sites[p] = displs_sites[p - 1] + recvcounts_sites[p - 1];
      numSitesGlobal += recvcounts_sites[p];
    }
  }

  // Global buffers (only on rank 0)
  std::vector<double> global_pos;
  std::vector<double> global_xi;
  std::vector<double> global_beta;
  std::vector<int> global_idx_beta_bcc;
  std::vector<int> global_rank;

  std::vector<int> recvcounts_pos, displs_pos;

  if (rank == 0) {
    global_pos.resize(dim * numSitesGlobal);
    global_xi.resize(numSitesGlobal);
    global_beta.resize(numSitesGlobal);
    global_idx_beta_bcc.resize(numSitesGlobal);
    global_rank.resize(numSitesGlobal);

    recvcounts_pos.resize(size);
    displs_pos.resize(size);
    for (int p = 0; p < size; ++p) {
      recvcounts_pos[p] = dim * recvcounts_sites[p];
      displs_pos[p] = dim * displs_sites[p];
    }
  }

  // Gatherv for positions
  MPI_Gatherv(local_pos.data(), dim * ns_local_int, MPI_DOUBLE,
              rank == 0 ? global_pos.data() : nullptr,
              rank == 0 ? recvcounts_pos.data() : nullptr,
              rank == 0 ? displs_pos.data() : nullptr, MPI_DOUBLE, 0, comm);

  // Helpers for scalar fields
  auto gatherv_scalar_double = [&](const std::vector<double>& local,
                                   std::vector<double>& global) {
    MPI_Gatherv(local.data(), ns_local_int, MPI_DOUBLE,
                rank == 0 ? global.data() : nullptr,
                rank == 0 ? recvcounts_sites.data() : nullptr,
                rank == 0 ? displs_sites.data() : nullptr, MPI_DOUBLE, 0, comm);
  };

  auto gatherv_scalar_int = [&](const std::vector<int>& local,
                                std::vector<int>& global) {
    MPI_Gatherv(local.data(), ns_local_int, MPI_INT,
                rank == 0 ? global.data() : nullptr,
                rank == 0 ? recvcounts_sites.data() : nullptr,
                rank == 0 ? displs_sites.data() : nullptr, MPI_INT, 0, comm);
  };

  gatherv_scalar_double(local_xi, global_xi);
  gatherv_scalar_double(local_beta, global_beta);
  gatherv_scalar_int(local_idx_beta_bcc, global_idx_beta_bcc);
  gatherv_scalar_int(local_rank, global_rank);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Rank 0 builds vtkPolyData and writes a single .vtp file
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (rank == 0) {
    auto points = vtkSmartPointer<vtkPoints>::New();

    auto atomic_number_arr = vtkSmartPointer<vtkDoubleArray>::New();
    atomic_number_arr->SetNumberOfComponents(1);
    atomic_number_arr->SetName("Atomic number");

    auto rho_arr = vtkSmartPointer<vtkDoubleArray>::New();
    rho_arr->SetNumberOfComponents(1);
    rho_arr->SetName("rho");

    auto beta_arr = vtkSmartPointer<vtkDoubleArray>::New();
    beta_arr->SetNumberOfComponents(1);
    beta_arr->SetName("thermalmultp");

    auto idx_beta_bcc_arr = vtkSmartPointer<vtkIntArray>::New();
    idx_beta_bcc_arr->SetNumberOfComponents(1);
    idx_beta_bcc_arr->SetName("thermalmultpbcc");

    auto rank_arr = vtkSmartPointer<vtkIntArray>::New();
    rank_arr->SetNumberOfComponents(1);
    rank_arr->SetName("rank");

    for (PetscInt i = 0; i < numSitesGlobal; ++i) {
      points->InsertNextPoint(global_pos[i * dim + 0], global_pos[i * dim + 1],
                              global_pos[i * dim + 2]);

      rho_arr->InsertNextValue(global_xi[i]);
      beta_arr->InsertNextValue(global_beta[i]);
      idx_beta_bcc_arr->InsertNextValue(global_idx_beta_bcc[i]);
      rank_arr->InsertNextValue(global_rank[i]);
    }

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(atomic_number_arr);
    polydata->GetPointData()->AddArray(rho_arr);
    polydata->GetPointData()->AddArray(beta_arr);
    polydata->GetPointData()->AddArray(idx_beta_bcc_arr);
    polydata->GetPointData()->AddArray(rank_arr);
    
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polydata);
#else
    writer->SetInputData(polydata);
#endif

    auto compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
    writer->SetCompressor(compressor);
    writer->Write();
  }

  MPI_Barrier(comm); 

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Restore all DMSwarm fields
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(
      DMSwarmRestoreField(atomistic_data, "idx", NULL, NULL, (void**)&idx_ptr));

  PetscCall(DMSwarmRestoreField(atomistic_data, "ghost", NULL, NULL,
                                (void**)&ghost_ptr));

  PetscCall(DMSwarmRestoreField(atomistic_data, DMSwarmPICField_coor, NULL,
                                NULL, (void**)&mean_q_ptr));

  PetscCall(DMSwarmRestoreField(atomistic_data, "molar-fraction", NULL, NULL,
                                (void**)&rho_ptr));

  PetscCall(DMSwarmRestoreField(atomistic_data, "beta", NULL, NULL,
                                (void**)&beta_ptr));

  PetscCall(DMSwarmRestoreField(atomistic_data, "idx-bcc-beta", NULL, NULL,
                                (void**)&idx_beta_bcc_ptr));

  PetscCall(DMSwarmRestoreField(atomistic_data, DMSwarmField_rank, NULL, NULL,
                                (void**)&rank_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/********************************************************************************/

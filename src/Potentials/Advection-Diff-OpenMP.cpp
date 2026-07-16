/**
 * @file MgHx-ADP-Residuals.cpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief
 * @version 0.1
 * @date 2025-10-14
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "petscsys.h"
#include "petscvec.h"
#include <cstdlib>
#if __APPLE__
#include <malloc/_malloc.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "Blobs/Blob-Function.hpp"
#include "Blobs/Blobs.hpp"
#include "Blobs/Topology.hpp"
#include "Macros.hpp"
#include "Mesh/Boundary-Conditions.hpp"
#include "petscvec.h"

extern PetscMPIInt size_MPI;
extern PetscMPIInt rank_MPI;

extern DiffusivePotential Potential_AD;

/********************************************************************************/

static PetscErrorCode evaluate_meassure(Vec rho,                             //!
                                        const Vec q_k1,                      //!
                                        const Vec beta_k1,                   //!
                                        const Vec mass,                      //!
                                        const ParticleTopology* atom_topology);  //!

static PetscErrorCode evaluate_F(PetscScalar* JKO_system,             //!
                                 PetscScalar Delta_t,                 //!
                                 const Vec rho,                       //!
                                 const Vec q_k1,                      //!
                                 const Vec q_k,                       //!
                                 const Vec beta_k1,                   //!
                                 const Vec beta_k,                    //!
                                 const Vec mass,                      //!
                                 const ParticleTopology* atom_topology);  //!

static PetscErrorCode evaluate_D_F_Dq(Vec D_JKO_Dq,                        //!
                                      PetscScalar Delta_t,                 //!
                                      const Vec rho_k1,                    //!
                                      const Vec x_k1,                      //!
                                      const Vec x_k,                       //!
                                      const Vec beta_k1,                   //!
                                      const Vec mass,                      //!
                                      const ParticleTopology* atom_topology);  //!

/********************************************************************************/

governing_equations Advection_Diff_constructor() {

  governing_equations equations;

  equations.evaluate_meassure_JKO = evaluate_meassure;

  equations.evaluate_JKO = evaluate_F;

  equations.evaluate_D_JKO_Dq = evaluate_D_F_Dq;

  return equations;
}

/************************************************************************/

static PetscErrorCode evaluate_meassure(
    Vec rho_k1,                         //!
    const Vec x_k1,                     //!
    const Vec beta_k1,                  //!
    const Vec mass,                     //!
    const ParticleTopology* atom_topology)  //!                        //!
{
  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  Blob shapefunc;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get local size
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof;
  PetscCall(VecGetLocalSize(x_k1, &n_dof));
  PetscInt n_sites = n_dof / dim;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to the local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec rho_k1_loc, mass_loc, x_k1_loc, x_k_loc, beta_k1_loc, beta_k_loc;
  PetscCall(VecGhostGetLocalForm(rho_k1, &rho_k1_loc));
  PetscCall(VecGhostGetLocalForm(mass, &mass_loc));
  PetscCall(VecGhostGetLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostGetLocalForm(beta_k1, &beta_k1_loc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to raw data from PETSc local (ghosted) vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar* rho_k1_ptr;
  PetscCall(VecGetArray(rho_k1_loc, &rho_k1_ptr));

  const PetscScalar* x_k1_ptr;
  PetscCall(VecGetArrayRead(x_k1_loc, &x_k1_ptr));

  const PetscScalar* beta_k1_ptr;
  PetscCall(VecGetArrayRead(beta_k1_loc, &beta_k1_ptr));

  const PetscScalar* mass_ptr;
  PetscCall(VecGetArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Evaluate the density for the Jordan-Kinderlerhrer-Otto functional
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_u = 0; site_u < n_sites; site_u++) {

    //! @brief Get atomistic information of site i
    PetscScalar rho_k1_u = 0.0;
    Eigen::Map<const Eigen::Vector3d> x_u_k1(&x_k1_ptr[dim * site_u]);

    //! Get topology data of subset Bu(rc)
    unsigned int size_Bu = atom_topology[site_u].size;
    const PetscInt* list_Bu = atom_topology[site_u].list;

    //! @brief Compute local sum
    for (int idx_i = 0; idx_i < size_Bu; idx_i++) {

      PetscInt site_i = list_Bu[idx_i];
      PetscScalar beta_i = beta_k1_ptr[site_i];
      PetscScalar m_i = mass_ptr[site_i];
      Eigen::Map<const Eigen::Vector3d> x_i_k1(&x_k1_ptr[dim * site_i]);

      rho_k1_u += shapefunc.N_i(x_u_k1, x_i_k1, beta_i) * m_i;
    }

    //! @brief Update local contribution of the residual equation
    rho_k1_ptr[site_u] = rho_k1_u;
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore raw data from PETSc local vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArray(rho_k1_loc, &rho_k1_ptr));
  PetscCall(VecRestoreArrayRead(x_k1_loc, &x_k1_ptr));
  PetscCall(VecRestoreArrayRead(beta_k1_loc, &beta_k1_ptr));
  PetscCall(VecRestoreArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostRestoreLocalForm(rho_k1, &rho_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(beta_k1, &beta_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(mass, &mass_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

static PetscErrorCode evaluate_F(
    PetscScalar* JKO_system,            //!
    PetscScalar Delta_t,                //!
    const Vec rho_k1,                   //!
    const Vec x_k1,                     //!
    const Vec x_k,                      //!
    const Vec beta_k1,                  //!
    const Vec beta_k,                   //!
    const Vec mass,                     //!
    const ParticleTopology* atom_topology)  //!                        //!
{
  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;
  PetscScalar kappa = Potential_AD.kappa;      //!
  PetscScalar rho_ref = Potential_AD.rho_ref;  //!

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get local size
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof;
  PetscCall(VecGetLocalSize(x_k1, &n_dof));
  PetscInt n_sites = n_dof / dim;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to the local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec rho_k1_loc, x_k1_loc, x_k_loc, beta_k1_loc, beta_k_loc, mass_loc;
  PetscCall(VecGhostGetLocalForm(rho_k1, &rho_k1_loc));
  PetscCall(VecGhostGetLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostGetLocalForm(x_k, &x_k_loc));
  PetscCall(VecGhostGetLocalForm(beta_k1, &beta_k1_loc));
  PetscCall(VecGhostGetLocalForm(beta_k, &beta_k_loc));
  PetscCall(VecGhostGetLocalForm(mass, &mass_loc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to raw data from PETSc local (ghosted) vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const PetscScalar* rho_k1_ptr;
  PetscCall(VecGetArrayRead(rho_k1_loc, &rho_k1_ptr));

  const PetscScalar* x_k1_ptr;
  PetscCall(VecGetArrayRead(x_k1_loc, &x_k1_ptr));

  const PetscScalar* x_k_ptr;
  PetscCall(VecGetArrayRead(x_k_loc, &x_k_ptr));

  const PetscScalar* beta_k1_ptr;
  PetscCall(VecGetArrayRead(beta_k1_loc, &beta_k1_ptr));

  const PetscScalar* beta_k_ptr;
  PetscCall(VecGetArrayRead(beta_k_loc, &beta_k_ptr));

  const PetscScalar* mass_ptr;
  PetscCall(VecGetArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Evaluate the Jordan-Kinderlerhrer-Otto functional
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  double JKO_local = 0.0;

#pragma omp parallel for reduction(+ : JKO_local) schedule(runtime)
  for (int site_i = 0; site_i < n_sites; site_i++) {

    //! @brief Get atomistic information of site i
    PetscScalar rho_k1_i = rho_k1_ptr[site_i];
    Eigen::Map<const Eigen::Vector3d> x_i_k1(&x_k1_ptr[dim * site_i]);
    Eigen::Map<const Eigen::Vector3d> x_i_k(&x_k_ptr[dim * site_i]);
    Eigen::Vector3d Dphi = x_i_k1 - x_i_k;
    PetscScalar m_i = mass_ptr[site_i];
    PetscScalar L2_Dphi = Dphi.norm();

    PetscScalar dr_i_k1 = 1.0 / sqrt(beta_k1_ptr[site_i]);
    PetscScalar dr_i_k = 1.0 / sqrt(beta_k_ptr[site_i]);
    PetscScalar Dr_i = dr_i_k1 - dr_i_k;

    //!
    PetscScalar JKO_i =
        (0.5 / Delta_t) * (L2_Dphi * L2_Dphi + Dr_i * Dr_i) * m_i +
        kappa * log(rho_k1_i / rho_ref) * m_i;

    //! @brief Update local contribution of the residual equation
    JKO_local += JKO_i;
  }

  //! Compute the resulting Meanfield Hamiltonian and Meanfield grand-cannonical
  //! partition function
  //! Perform partial sum reduction of each MPI process
  PetscCall(MPIU_Allreduce(&JKO_local, JKO_system, 1, MPI_DOUBLE, MPI_SUM,
                           MPI_COMM_WORLD));

  //! Check if the internal energy is a NAN
  if (PetscIsNanReal(*JKO_system) == PETSC_TRUE) {
    PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "evaluate_F", __FILE__,
                         PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                         "The internal energy takes a NaN value"));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  //! Check if the internal energy is a inft
  if (PetscIsInfReal(*JKO_system) == PETSC_TRUE) {
    PetscCall(PetscError(PETSC_COMM_SELF, __LINE__, "evaluate_F", __FILE__,
                         PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                         "The internal energy takes a infinity value"));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore raw data from PETSc local vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArrayRead(rho_k1_loc, &rho_k1_ptr));
  PetscCall(VecRestoreArrayRead(x_k1_loc, &x_k1_ptr));
  PetscCall(VecRestoreArrayRead(x_k_loc, &x_k_ptr));
  PetscCall(VecRestoreArrayRead(beta_k1_loc, &beta_k1_ptr));
  PetscCall(VecRestoreArrayRead(beta_k_loc, &beta_k_ptr));
  PetscCall(VecRestoreArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostRestoreLocalForm(rho_k1, &rho_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(x_k, &x_k_loc));
  PetscCall(VecGhostRestoreLocalForm(beta_k1, &beta_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(beta_k, &beta_k_loc));
  PetscCall(VecGhostRestoreLocalForm(mass, &mass_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

static PetscErrorCode evaluate_D_F_Dq(Vec D_JKO_Dq,                       //!
                                      PetscScalar Delta_t,                //!
                                      const Vec rho_k1,                   //!
                                      const Vec x_k1,                     //!
                                      const Vec x_k,                      //!
                                      const Vec beta_k1,                  //!
                                      const Vec mass,                     //!
                                      const ParticleTopology* atom_topology)  //!

{
  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  PetscScalar kappa = Potential_AD.kappa;      //!
  PetscScalar rho_ref = Potential_AD.rho_ref;  //!

  Blob shapefunc;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get local size
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof;
  PetscCall(VecGetLocalSize(x_k1, &n_dof));
  PetscInt n_sites = n_dof / dim;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to the local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec D_JKO_Dq_loc, rho_k1_loc, x_k1_loc, x_k_loc, beta_k1_loc, mass_loc;
  PetscCall(VecGhostGetLocalForm(D_JKO_Dq, &D_JKO_Dq_loc));
  PetscCall(VecGhostGetLocalForm(rho_k1, &rho_k1_loc));
  PetscCall(VecGhostGetLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostGetLocalForm(x_k, &x_k_loc));
  PetscCall(VecGhostGetLocalForm(beta_k1, &beta_k1_loc));
  PetscCall(VecGhostGetLocalForm(mass, &mass_loc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to raw data from PETSc local (ghosted) vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar* D_JKO_Dq_ptr;
  PetscCall(VecGetArray(D_JKO_Dq_loc, &D_JKO_Dq_ptr));

  const PetscScalar* rho_k1_ptr;
  PetscCall(VecGetArrayRead(rho_k1_loc, &rho_k1_ptr));

  const PetscScalar* x_k1_ptr;
  PetscCall(VecGetArrayRead(x_k1_loc, &x_k1_ptr));

  const PetscScalar* x_k_ptr;
  PetscCall(VecGetArrayRead(x_k_loc, &x_k_ptr));

  const PetscScalar* beta_k1_ptr;
  PetscCall(VecGetArrayRead(beta_k1_loc, &beta_k1_ptr));

  const PetscScalar* mass_ptr;
  PetscCall(VecGetArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Evaluate the Jordan-Kinderlerhrer-Otto functional
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_u = 0; site_u < n_sites; site_u++) {

    //! @brief Get atomistic information of site i
    PetscScalar rho_k1_u = rho_k1_ptr[site_u];
    PetscScalar m_u = mass_ptr[site_u];
    Eigen::Vector3d d_rho_k1_dx_u = Eigen::Vector3d::Zero();
    Eigen::Map<const Eigen::Vector3d> x_u_k1(&x_k1_ptr[dim * site_u]);
    Eigen::Map<const Eigen::Vector3d> x_u_k(&x_k_ptr[dim * site_u]);
    Eigen::Vector3d Dphi = x_u_k1 - x_u_k;

    //! Get topology data of subset Bu(rc)
    unsigned int size_Bu = atom_topology[site_u].size;
    const PetscInt* list_Bu = atom_topology[site_u].list;

    //! @brief Compute local gradient of the density
    for (int idx_i = 1; idx_i < size_Bu; idx_i++) {

      PetscInt site_i = list_Bu[idx_i];
      PetscScalar beta_i = beta_k1_ptr[site_i];
      PetscScalar m_i = mass_ptr[site_i];
      Eigen::Map<const Eigen::Vector3d> x_i_k1(&x_k1_ptr[dim * site_i]);

      d_rho_k1_dx_u += shapefunc.dN_i(x_u_k1, x_i_k1, beta_i) * m_i;
    }

    //!
    Eigen::Vector3d D_JKO_Dq_u =
        (Dphi / Delta_t) * m_u + kappa * (d_rho_k1_dx_u / rho_k1_u) * m_u;

    //! Fill residual vector
    for (PetscInt alpha = 0; alpha < dim; alpha++) {
      D_JKO_Dq_ptr[site_u * dim + alpha] = D_JKO_Dq_u(alpha);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore raw data from PETSc local vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArray(D_JKO_Dq_loc, &D_JKO_Dq_ptr));
  PetscCall(VecRestoreArrayRead(rho_k1_loc, &rho_k1_ptr));
  PetscCall(VecRestoreArrayRead(x_k1_loc, &x_k1_ptr));
  PetscCall(VecRestoreArrayRead(x_k_loc, &x_k_ptr));
  PetscCall(VecRestoreArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostRestoreLocalForm(D_JKO_Dq, &D_JKO_Dq_loc));
  PetscCall(VecGhostRestoreLocalForm(rho_k1, &rho_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(x_k, &x_k_loc));
  PetscCall(VecGhostRestoreLocalForm(beta_k1, &beta_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(mass, &mass_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/
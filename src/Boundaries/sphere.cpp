

#include "Boundaries/sphere.hpp"

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
#include "Macros.hpp"

extern PetscMPIInt size_MPI;
extern PetscMPIInt rank_MPI;

/************************************************************************/

PetscErrorCode BC_Sphere::add_barrier_potential(PetscScalar *JKO_system,
                                                const Vec x_k1,
                                                const Vec mass)
{
  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get local size
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof;
  PetscCall(VecGetLocalSize(x_k1, &n_dof));
  PetscInt n_sites = n_dof / dim;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to the local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec x_k1_loc, mass_loc;
  PetscCall(VecGhostGetLocalForm(mass, &mass_loc));
  PetscCall(VecGhostGetLocalForm(x_k1, &x_k1_loc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to raw data from PETSc local (ghosted) vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  const PetscScalar *x_k1_ptr;
  PetscCall(VecGetArrayRead(x_k1_loc, &x_k1_ptr));

  const PetscScalar *mass_ptr;
  PetscCall(VecGetArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Evaluate the barrier potential
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  double barrier_potential_local = 0.0;
  double barrier_potential = 0.0;

#pragma omp parallel for reduction(+ : barrier_potential_local) \
    schedule(runtime)
  for (int site_u = 0; site_u < n_sites; site_u++)
  {

    //! @brief Get particle information of site i
    Eigen::Map<const Eigen::Vector3d> x_u_k1(&x_k1_ptr[dim * site_u]);
    PetscScalar mass_u = mass_ptr[site_u];

    //! @brief Compute the distance between the particle and the center of the
    //! sphere
    PetscScalar distance = (x_u_k1 - center_sphere).norm() - radius_sphere;

    //! @brief Add the barrier potential to the JKO system
    barrier_potential_local += 0.5 * penalty * DSQR(distance) * mass_u;
  }

  //! Perform partial sum reduction of each MPI process
  PetscCall(MPIU_Allreduce(&barrier_potential_local, &barrier_potential, 1,
                           MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

  //! Check if the barrier potential is a NAN
  if (PetscIsNanReal(barrier_potential) == PETSC_TRUE)
  {
    PetscCall(PetscError(PETSC_COMM_SELF, __LINE__,
                         "BC_Sphere::add_barrier_potential", __FILE__,
                         PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                         "The barrier potential takes a NaN value"));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  //! Check if the barrier potential is infinity
  if (PetscIsInfReal(barrier_potential) == PETSC_TRUE)
  {
    PetscCall(PetscError(PETSC_COMM_SELF, __LINE__,
                         "BC_Sphere::add_barrier_potential", __FILE__,
                         PETSC_ERR_RETURN, PETSC_ERROR_INITIAL,
                         "The barrier potential takes a infinity value"));
    PetscFunctionReturn(PETSC_ERR_RETURN);
  }

  //! Update the barrier potential
  *JKO_system += barrier_potential;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore raw data from PETSc local vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArrayRead(x_k1_loc, &x_k1_ptr));
  PetscCall(VecRestoreArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostRestoreLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(mass, &mass_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/

PetscErrorCode BC_Sphere::add_barrier_forces(Vec D_JKO_Dq, const Vec x_k1,
                                             const Vec mass)
{

  PetscFunctionBeginUser;

  unsigned int dim = NumberDimensions;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Get local size
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt n_dof;
  PetscCall(VecGetLocalSize(x_k1, &n_dof));
  PetscInt n_sites = n_dof / dim;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to the local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  Vec D_JKO_Dq_loc, mass_loc, x_k1_loc;
  PetscCall(VecGhostGetLocalForm(D_JKO_Dq, &D_JKO_Dq_loc));
  PetscCall(VecGhostGetLocalForm(mass, &mass_loc));
  PetscCall(VecGhostGetLocalForm(x_k1, &x_k1_loc));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Acces to raw data from PETSc local (ghosted) vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar *D_JKO_Dq_ptr;
  PetscCall(VecGetArray(D_JKO_Dq_loc, &D_JKO_Dq_ptr));

  const PetscScalar *x_k1_ptr;
  PetscCall(VecGetArrayRead(x_k1_loc, &x_k1_ptr));

  const PetscScalar *mass_ptr;
  PetscCall(VecGetArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Evaluate the barrier forces
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#pragma omp parallel for schedule(runtime)
  for (int site_u = 0; site_u < n_sites; site_u++)
  {

    //! @brief Get particle information of site i
    Eigen::Map<const Eigen::Vector3d> x_u_k1(&x_k1_ptr[dim * site_u]);
    PetscScalar mass_u = mass_ptr[site_u];

    //! @brief Compute the distance between the particle and the center of the
    //! sphere
    PetscScalar distance = (x_u_k1 - center_sphere).norm() - radius_sphere;

    //! @brief Compute the barrier forceforce_u
    Eigen::Vector3d force_u = - penalty * distance * (x_u_k1 - center_sphere).normalized() * mass_u;
    //! Fill residual vector
    for (PetscInt alpha = 0; alpha < dim; alpha++)
    {
      D_JKO_Dq_ptr[site_u * dim + alpha] += force_u(alpha);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore raw data from PETSc local vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecRestoreArray(D_JKO_Dq_loc, &D_JKO_Dq_ptr));
  PetscCall(VecRestoreArrayRead(x_k1_loc, &x_k1_ptr));
  PetscCall(VecRestoreArrayRead(mass_loc, &mass_ptr));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Restore local version of the vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGhostRestoreLocalForm(D_JKO_Dq, &D_JKO_Dq_loc));
  PetscCall(VecGhostRestoreLocalForm(x_k1, &x_k1_loc));
  PetscCall(VecGhostRestoreLocalForm(mass, &mass_loc));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/************************************************************************/
/**
 * @file Macros.hpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief
 * @version 0.1
 * @date 2022-07-19
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _MACROS_H_
#define _MACROS_H_

#include "petscis.h"
#include <Eigen/Dense>
#include <iostream>
#include <petsc/private/dmimpl.h>
#include <petscao.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmshell.h>
#include <petscdmswarm.h>
#include <petscerror.h>
#include <petscsf.h>
#include <petscsystypes.h>
#include <petscvec.h>

#ifndef PETSC_SUCCESS
#define PETSC_SUCCESS 0
#endif

#ifndef PETSC_ERR_NOT_CONVERGED
#define PETSC_ERR_NOT_CONVERGED 91 // solver did not converge
#endif

#ifndef PETSC_ERR_RETURN
#define PETSC_ERR_RETURN 99
#endif

#define sqrt_2 1.4142135623730951
#define sqrt_3 1.7320508075688772
#define sqrt_6 2.449489742783178

//! Auxiliar Eigen function to map
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> VectorType;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    MatrixType;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> List1D;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    List2D;

/*******************************************************/

/**
 * @brief  Ising model type
 *
 */
enum class BackgroundMeshType {
  DMDA_mesh,  //!< DMDA mesh
  DMPLEX_mesh //!< DMPLEX mesh
};

/*******************************************************/

enum Miller_Index {
  //! @param idx_0_0_0: direction (0, 0, 0)
  idx_0_0_0,
  //! @param idx_m1_m1_m1: direction (-1, -1, -1)
  idx_m1_m1_m1,
  //! @param idx_m1_m1_0: direction (-1, -1, 0)
  idx_m1_m1_0,
  //! @param idx_m1_m1_0: direction (-1, -1, 1)
  idx_m1_m1_p1,
  //! @param idx_m1_0_m1: direction (-1, 0, -1)
  idx_m1_0_m1,
  //! @param idx_m1_0_0: direction (-1, 0, 0)
  idx_m1_0_0,
  //! @param idx_m1_0_p1: direction (-1, 0, 1)
  idx_m1_0_p1,
  //! @param idx_m1_p1_m1: direction (-1, 1, -1)
  idx_m1_p1_m1,
  //! @param idx_m1_p1_0: direction (-1, 1, 0)
  idx_m1_p1_0,
  //! @param idx_m1_p1_p1: direction (-1, 1, 1)
  idx_m1_p1_p1,
  //! @param idx_0_m1_m1: direction (0, -1, -1)
  idx_0_m1_m1,
  //! @param idx_0_m1_0: direction (0, -1, 0)
  idx_0_m1_0,
  //! @param idx_0_m1_p1: direction (0, -1, 1)
  idx_0_m1_p1,
  //! @param idx_0_0_m1: direction (0, 0, -1)
  idx_0_0_m1,
  //! @param idx_0_0_p1: direction (0, 0, 1)
  idx_0_0_p1,
  //! @param idx_0_p1_m1: direction (0, 1, -1)
  idx_0_p1_m1,
  //! @param idx_0_p1_0: direction (0, 1, 0)
  idx_0_p1_0,
  //! @param idx_0_p1_p1: direction (0, 1, 1)
  idx_0_p1_p1,
  //! @param idx_p1_m1_m1: direction (1, -1, -1)
  idx_p1_m1_m1,
  //! @param idx_p1_m1_0: direction (1, -1, 0)
  idx_p1_m1_0,
  //! @param idx_p1_m1_p1: direction (1, -1, 1)
  idx_p1_m1_p1,
  //! @param idx_p1_0_m1: direction (1, 0, -1)
  idx_p1_0_m1,
  //! @param idx_p1_0_0: direction (1, 0, 0)
  idx_p1_0_0,
  //! @param idx_p1_0_p1: direction (1, 0, 1)
  idx_p1_0_p1,
  //! @param idx_p1_p1_m1: direction (1, 1, -1)
  idx_p1_p1_m1,
  //! @param idx_p1_p1_0: direction (1, 1, 0)
  idx_p1_p1_0,
  //! @param idx_p1_p1_p1: direction (1, 1, 1)
  idx_p1_p1_p1,
};

/*******************************************************/

typedef struct {

  PetscInt size;

  const PetscInt *list;

} ParticleTopology;

/*******************************************************/

typedef struct dump_file {

  /** @param n_particles Number of particles */
  int n_particles;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Box bounds
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  double box_x_min, box_x_max;
  double box_y_min, box_y_max;
  double box_z_min, box_z_max;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Boundary conditions
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMBoundaryType bx;
  DMBoundaryType by;
  DMBoundaryType bz;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    System Lagrange multipliers
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /** @param beta: Thermodynamic Lagrange multiplier */
  double *beta;

  /** @param beta_bcc: Index for the sites with beta boundary condition */
  int *beta_bcc;

  /** @param gamma: Chemical Lagrange multiplier */
  double *gamma;

  /** @param gamma_bcc: Index for the sites with gamma boundary condition */
  int *gamma_bcc;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Position-related variables
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /** @param mean_q: Mean value of each particle position */
  double *mean_q;

  /** @param stdv_q: Standard desviation of exach particle position */
  double *stdv_q;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Chemical variables
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /** @param xi: Molar fraction (mean occupancy) */
  double *xi;

} dump_file;

/*******************************************************/

/**
 * @brief Abstract base class for user-defined JKO / mass-transport equations.
 *
 * Derive from this class to provide a concrete potential (e.g. advection–
 * diffusion). Call sites should pass a reference (or non-owning pointer) so
 * that polymorphic overrides are invoked.
 */
class GoverningEquations {
public:
  virtual ~GoverningEquations() {}

  /**
   * @brief Reconstruct the discrete density measure from particle data.
   */
  virtual PetscErrorCode
  evaluate_meassure_JKO(Vec rho, const Vec q_k1, const Vec beta_k1,
                        const Vec mass,
                        const ParticleTopology *particle_topology) = 0;

  /**
   * @brief Evaluate the discrete Jordan–Kinderlehrer–Otto functional
   *
   \f[
    F(\rho_{k+1}) \sim \sum_p^n \frac{1}{2} \frac{\|q_{p,k+1} - q_{p,k} \|^2 +
   \|\beta_{p,k+1} - \beta_{p,k} \|^2}{t_{k+1} - t_{k}} m_p +
   V(\rho_{k+1}(x_{p,k+1}))\, m_p
   \f]
   */
  virtual PetscErrorCode
  evaluate_JKO(PetscScalar *JKO_system, PetscScalar Delta_t, const Vec rho,
               const Vec q_k1, const Vec q_k, const Vec beta_k1,
               const Vec beta_k, const Vec mass,
               const ParticleTopology *particle_topology) = 0;

  /**
   * @brief Evaluate the position gradient of the discrete JKO functional.
   */
  virtual PetscErrorCode
  evaluate_D_JKO_Dq(Vec D_JKO_Dq, PetscScalar Delta_t, const Vec rho_k1,
                    const Vec x_k1, const Vec x_k, const Vec beta_k1,
                    const Vec mass,
                    const ParticleTopology *particle_topology) = 0;
};

/*******************************************************/

class boundaryCondition {
public:
  virtual ~boundaryCondition() {}

  /**
   * @brief Add the barrier potential
   */
  virtual PetscErrorCode add_barrier_potential(PetscScalar *JKO_system,
                                               const Vec x_k1,
                                               const Vec mass) = 0;

  /**
   * @brief Add the barrier forces to the JKO system.
   */
  virtual PetscErrorCode add_barrier_forces(Vec D_JKO_Dq, const Vec x_k1,
                                            const Vec mass) = 0;
};

/*******************************************************/

/**
 * @brief Abstract base class for Shape Functions and Potential Components.
 * * This interface defines the contract for any function N_i and its
 * spatial gradient DN_i, evaluated at point x relative to a reference
 * position x_i (e.g., a node or an particle site).
 */
class ShapeFunction {
public:
  /**
   * @brief Virtual destructor to ensure proper cleanup of derived classes.
   */
  virtual ~ShapeFunction() {}

  /**
   * @brief Evaluates the function value at position x relative to x_i.
   * * @param x The evaluation point coordinates.
   * @param x_i The reference position or nodal coordinates.
   * @return The scalar value of the function.
   */
  virtual double N_i(const Eigen::Vector3d &x, const Eigen::Vector3d &x_i,
                     double beta_i) const = 0;

  /**
   * @brief Calculates the gradient (first derivative) of the function.
   * * @param x The evaluation point coordinates.
   * @param x_i The reference position or nodal coordinates.
   * @return A Vector3d containing the partial derivatives [dN/dx, dN/dy,
   * dN/dz].
   */
  virtual Eigen::Vector3d dN_i(const Eigen::Vector3d &x,
                               const Eigen::Vector3d &x_i,
                               double beta_i) const = 0;
};

/*******************************************************/

/*
  Color text
*/
#define RESET "\033[0m"
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

/*
  SOLERA units and constants
*/
// Reduced Planck constant -> 6.5821192815×10−4 [eV·ps]
#define h_planck 6.5821192815E-4
// Boltzmann constant [ev/K]
#define k_B 8.617332478E-5
// 1 u.m.a = 103.4993333E-6 [ev·A^{-2}·ps^{2}]
#define unit_change_uma 103.4993333E-6
// the factor 98.22694969 is to change the
// units: 98.22694969 1/ps = 1 ((eV/A^2)/u.m.a.)^(1/2)
#define unit_change_w 98.227002603

/*
  Constant macros
*/
#define MAXW 100
#define MAXC 1000
#define NumberDimensions 3
#define TOL_NR 10E-6
#define TOL_zero 10E-23
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//! The length of the buffer used to efficient dynamic re-sizing
#define BufferLenght 100

//! Angular dependant potential cutoff
#define r_cutoff_default 6

#define maxneigh 500

/*
  Math macros
*/
static float sqr_arg;
#define SQR(a) ((sqr_arg = (a)) == 0.0 ? 0.0 : sqr_arg * sqr_arg)
static double dsqr_arg;
#define DSQR(a) ((dsqr_arg = (a)) == 0.0 ? 0.0 : dsqr_arg * dsqr_arg)
static double dmax_arg1, dmax_arg2;
#define DMAX(a, b)                                                             \
  (dmax_arg1 = (a), dmax_arg2 = (b),                                           \
   (dmax_arg1) > (dmax_arg2) ? (dmax_arg1) : (dmax_arg2))
static double dmin_arg1, dmin_arg2;
#define DMIN(a, b)                                                             \
  (dmin_arg1 = (a), dmin_arg2 = (b),                                           \
   (dmin_arg1) < (dmin_arg2) ? (dmin_arg1) : (dmin_arg2))
static float max_arg1, max_arg2;
#define FMAX(a, b)                                                             \
  (max_arg1 = (a), max_arg2 = (b),                                             \
   (max_arg1) > (max_arg2) ? (max_arg1) : (max_arg2))
static float min_arg1, min_arg2;
#define FMIN(a, b)                                                             \
  (min_arg1 = (a), min_arg2 = (b),                                             \
   (min_arg1) < (min_arg2) ? (min_arg1) : (min_arg2))
static long lmax_arg1, lmax_arg2;
#define LMAX(a, b)                                                             \
  (lmax_arg1 = (a), lmax_arg2 = (b),                                           \
   (lmax_arg1) > (lmax_arg2) ? (lmax_arg1) : (lmax_arg2))
static long lmin_arg1, lmin_arg2;
#define LMIN(a, b)                                                             \
  (lmin_arg1 = (a), lmin_arg2 = (b),                                           \
   (lmin_arg1) < (lmin_arg2) ? (lmin_arg1) : (lmin_arg2))
static int imax_arg1, imax_arg2;
#define IMAX(a, b)                                                             \
  (imax_arg1 = (a), imax_arg2 = (b),                                           \
   (imax_arg1) > (imax_arg2) ? (imax_arg1) : (imax_arg2))
static int imin_arg1, imin_arg2;
#define IMIN(a, b)                                                             \
  (imin_arg1 = (a), imin_arg2 = (b),                                           \
   (imin_arg1) < (imin_arg2) ? (imin_arg1) : (imin_arg2))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) s

#endif
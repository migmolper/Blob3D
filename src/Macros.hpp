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
#define PETSC_ERR_NOT_CONVERGED 91  // solver did not converge
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

/**
 * @brief Create dictionary for the species
 *
 */
enum AtomicSpecie {
  //! @param Rigid: Rigid body
  Rigid,
  //! @param H: Hydrogen, Z: 1.0079
  H,
  //! @param He: Helium, Z: 4.0026
  He,
  Li,
  Be,
  B,
  C,
  N,
  O,
  F,
  Ne,
  Na,
  //! @param Mg: Magnesium, Z: 24.305
  Mg,
  //! @param Al: Aluminium, Z: 26.981
  Al,
  Si,
  P,
  S,
  Cl,
  Ar,
  K,
  Ca,
  Sc,
  Ti,
  V,
  Cr,
  Mn,
  Fe,
  Co,
  Ni,
  //! @param Cu: Copper, Z: 63,546
  Cu,
  Zn,
  Ga,
  Ge,
  As,
  Se,
  Br,
  Kr,
  Rb,
  Sr,
  Y,
  Zr,
  Nb,
  Mo,
  Tc,
  Ru,
  Rh,
  Pd,
  Ag,
  Cd,
  In,
  Sn,
  Sb,
  Te,
  I,
  Xe,
  Cs,
  Ba,
  La,
  Ce,
  Pr,
  Nd,
  Pm,
  Sm,
  Eu,
  Gd,
  Tb,
  Dy,
  Ho,
  Er,
  Tm,
  Yb,
  Lu,
  Hf,
  Ta,
  //! @param W: Wolfram (Tungsten), Z: 183.84
  W,
  Re,
  Os,
  Ir,
  Pt,
  Au,
  Hg,
  Tl,
  Pb,
  Bi,
  Po,
  At,
  Rn,
  Fr,
  Ra,
  Ac,
  Th,
  Pa,
  U,
  Np,
  Pu,
  Am,
  Cm,
  Bk,
  Cf,
  Es,
  Fm,
  Md,
  No,
  Lr,
  Rf,
  Db,
  Sg,
  Bh,
  Hs,
  Mt,
  Ds,
  Rg
};

/*******************************************************/

enum species_combinations {
  HH,
  AlAl,
  CuCu,
  FeFe,
  MgMg,
  FeH,
  HFe,
  MgH,
  HMg,
  AlCu,
  CuAl
};

/*******************************************************/

/**
 * @brief  Ising model type
 *
 */
enum class BackgroundMeshType {
  DMDA_mesh,   //!< DMDA mesh
  DMPLEX_mesh  //!< DMPLEX mesh
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

  const PetscInt* list;

} AtomTopology;

/*******************************************************/

/**
 * @brief This structures defines a function and its derivatives
 *
 */
typedef struct {

  /*! @param F: Integrand */
  void (*F)(double* F, const double* xi, const double* q,
            const AtomicSpecie* spc);

  /*! @param dF_dq: Gradient of the function (analytical) */
  void (*dF_dq)(int direction, double* dF_dq, const double* xi, const double* q,
                const AtomicSpecie* spc);

  /*! @param d2F_dq2: Hessian of the function (analytical) */
  void (*d2F_dq2)(int direction, double* d2F_dq2, const double* xi,
                  const double* q, const AtomicSpecie* spc);

  /*! @param dF_dq_FD: Gradient of the function (numerical) */
  void (*dF_dq_FD)(int direction, double* dF_dq, const double* xi,
                   const double* q, const AtomicSpecie* spc);

  /*! @param d2F_dq2_FD: Hessian of the function (numerical) */
  void (*d2F_dq2_FD)(int direction, double* d2F_dq2, const double* xi,
                     const double* q, const AtomicSpecie* spc);

  /*! @param dF_dn: Gradient of the function with respect the occupancy */
  void (*dF_dn)(int direction, double* dF_dn, const double* xi, const double* q,
                const AtomicSpecie* spc);

} potential_function;

/*******************************************************/

/**
 * @brief Nonuniform cubic splines with n intervals
 *
 */
typedef struct CubicSpline {
  double dx;
  int n;
  double* x;
  double* a;
  double* b;
  double* c;
  double* d;
  double* db;
  double* dc;
  double* dd;
  double* ddc;
  double* ddd;
} CubicSpline;

/*******************************************************/

typedef struct dump_file {

  /** @param n_atoms Number of atoms */
  int n_atoms;

  /** @param specie: Integer which defines the atomic specie
   */
  AtomicSpecie* specie;

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
  double* beta;

  /** @param beta_bcc: Index for the sites with beta boundary condition */
  int* beta_bcc;

  /** @param gamma: Chemical Lagrange multiplier */
  double* gamma;

  /** @param gamma_bcc: Index for the sites with gamma boundary condition */
  int* gamma_bcc;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Position-related variables
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /** @param mean_q: Mean value of each atomic position */
  double* mean_q;

  /** @param stdv_q: Standard desviation of each atomic position */
  double* stdv_q;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Chemical variables
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /** @param xi: Molar fraction (mean occupancy) */
  double* xi;

} dump_file;

/*******************************************************/

/**
 * @brief  Ising model type
 *
 */
enum class IsingModelType {
  LatticeGas,  //!< Lattice gas model (interestitial)
  BinaryAlloy  //!< Binary alloy model (substitutional)
};

/*******************************************************/

/**
 * @brief Structure which contains the necessary information to evaluate the
 * ADP interatomic potential
 *
 */
typedef struct adpPotential {

  /**
   * @brief
   *
   */
  int n_embed;

  int n_rho;
  int n_pair;
  int n_u;
  int n_w;

  /**
   * @brief Evaluation of the embedded energy
   *
   */
  CubicSpline embed;

  /**
   * @brief Evaluation of the density function
   *
   */
  CubicSpline rho;

  /**
   * @brief
   *
   */
  CubicSpline pair;
  CubicSpline u;
  CubicSpline w;

  /**
   * @brief Mass of the specie
   *
   */
  double mass;
  double radius;
  double factor;

  /**
   * @brief Cut-off radious for the specie
   *
   */
  double r_cutoff;

} adpPotential;

/********************************************************************************/

/**
 * @brief Structure which contains the necessary information to evaluate the
 * EAM interatomic potential
 *
 */
typedef struct eamPotential {

  /**
   * @brief
   *
   */
  int n_embed;

  int n_rho;
  int n_pair;

  /**
   * @brief Evaluation of the embedded energy
   *
   */
  CubicSpline embed;

  /**
   * @brief Evaluation of the density function
   *
   */
  CubicSpline rho;

  /**
   * @brief
   *
   */
  CubicSpline pair;

  /**
   * @brief Mass of the specie
   *
   */
  double mass;
  double radius;
  double factor;

  /**
   * @brief Cut-off radious for the specie
   *
   */
  double r_cutoff;

} eamPotential;

/********************************************************************************/

/**
 * @brief Class which contains the necessary information to evaluate the
 * indenter potential.
 * \f[
 *  U = \varepsilon \left( \frac{\sigma}{r} \right)^{12}
 * \f]
 */
class indenterPotential {

 private:
  double sigma;     //!< Repulsive potential strength
  double epsilon;   //!< Repulsive potential depth
  double r_cutoff;  //!< Repulsive potential cutoff

  void dV_indenter_ij_dq(double* dV_dq, const double* n, const double* q,
                         const AtomicSpecie* spc) const;

 public:
  indenterPotential();

  /**
   * @brief Initialize the indenter potential parameters.
   * @param sigma: Repulsive potential strength (nullptr uses default).
   * @param epsilon: Repulsive potential depth (nullptr uses default).
   * @param r_cutoff: Repulsive potential cutoff (nullptr uses default).
   */
  PetscErrorCode init(const double* sigma = nullptr,
                      const double* epsilon = nullptr,
                      const double* r_cutoff = nullptr);

  /**
   * @brief Evaluate the gradient of the indenter potential with respect to the
   * mean atomic positions
   * @param RHS: Residual / force vector (mean-q block), updated in place.
   * @param mean_q: Mean atomic positions (ghosted).
   * @param xi: Site occupancies (ghosted).
   * @param specie_ptr: Local species array (ghosted layout).
   * @param atom_topology: Neighbour lists.
   * @return PetscErrorCode PETSC_SUCCESS if successful, PETSC_ERR_OTHER if
   * error.
   */
  PetscErrorCode evaluate_dV_indenter_dmq(Vec RHS, const Vec mean_q,
                                          const Vec xi,
                                          const AtomicSpecie* specie_ptr,
                                          const AtomTopology* atom_topology) const;
};

/********************************************************************************/


/**
 * @brief Global variable to define a DMD simulation
 *
 */
 typedef struct DMD {

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   @brief System information
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  //! @param Variable: with the atomistic data
  DM atomistic_data;

  //! @param IsingModel: Ising model (interstitial or subtitutional)
  IsingModelType IsingModel;

  //! @param n_sites_global: Number of atoms in the global domain
  PetscInt n_sites_global;

  //! @param n_sites_local: Number of atoms in the local domain (without ghost)
  PetscInt n_sites_local;

  //! @param n_ghost: Number of ghost atoms in the local domain
  PetscInt n_ghost;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  @brief Enviroment thermo-chemo-mechanical variables
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  //! @param Pressure_env: Enviromental value of the pressure
  double Pressure_env;

  //! @param Temperature_env: Enviromental value of the temperature
  double Temperature_env;

  //! @param ChemicalPotential_env: Enviromental value of the chemical potential
  //! ({mu}). {mu} = {gamma}/{beta}
  double ChemicalPotential_env;

  //! @param MolarFraction_env: Enviromental value of the chemical potential
  //! ({mu}). {mu} = {gamma}/{beta}
  double MolarFraction_env;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @brief Transition state variables
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  //! @param Interstitial
  double particle_mass_I;

  //! @param Subtitutional
  double particle_mass_A;
  double particle_mass_B;

  //! @param
  double max_diffusion_lenght;

  //! @param
  double Eb;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   @brief Topological variables of each site
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  /** @param neigh: Local indices of I atoms (we only need for the cell) */
  IS* mechanical_neighs_idx;

  /** @param diffusive_neighs_idx: Table with the list of neighbors-idx of
   * each diffusive site */
  IS* diffusive_neighs_idx;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @brief Topological variables for the themo-mechanical equations
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  /** @param n_mechanical_sites_local: Number of mechanical sites (local) */
  PetscInt n_mechanical_sites_local;

  /** @param n_mechanical_sites_ghost: Number of mechanical ghost sites */
  PetscInt n_mechanical_sites_ghost;

  /** @param active_mech_sites: List with the active mechanical sites */
  IS active_mech_sites;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @brief Topological variables for the chemical equation
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  /** @param n_active_diff_sites_local: Number of active diffusive sites (local)
   */
  PetscInt n_active_diff_sites_local;

  /** @param n_active_diff_sites_ghost: Number of active diffusive ghost sites
   */
  PetscInt n_active_diff_sites_ghost;

  /** @param active_diff_sites: List with the active diffusive sites using local
   * numbering */
  IS active_diff_sites;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @brief Miscelaneous variables
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  /** @param dump2petsc_mapping: Allow the user to track the index of each atom
   * from the .dump prdering and petsc ordering */
  AO dump2petsc_mapping;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @brief Indenter potential
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  //! @param V_indenter: Indenter potential
  indenterPotential V_indenter;

} DMD;

/*******************************************************/

/*******************************************************/

/**
 * @brief Structure which contains the necessary information to evaluate the
 * DiffusivePotential potential
 *
 */
typedef struct DiffusivePotential {

  //! @brief Diffusion coefficient
  PetscScalar kappa;

  //! @brief Reference density of the system
  PetscScalar rho_ref;

  //! @brief Penaly stiffness
  PetscScalar C;

} DiffusivePotential;

/*******************************************************/

/**
 * @brief Structure which contains user define equations to evaluate any sort of
 * equilibrium equation
 *
 */
typedef struct dmd_equations {

  /**
   * @brief Evaluate the local density of the system
   *
   * @param mean_q:Mean value of q
   * @param xi Molar fraction
   * @param specie Atomic specie
   * @param atom_topology_i List of neighs
   * @return rho
   */
  PetscErrorCode (*evaluate_rho)(Vec rho,                             //!
                                 const Vec mean_q,                    //!
                                 const Vec xi,                        //!
                                 const AtomicSpecie* specie_ptr,      //!
                                 const AtomTopology* atom_topology);  //!

  /**
   * @brief Evaluate the local potential of the system
   * @param mean_q Mean position (\f$\bar{\mathbf{q}}\f$) of the atomic
   * positions in the updated configuration
   * @param xi Molar fraction
   * @param rho Local density of the system
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return V_system Total potential of the system
   */
  PetscErrorCode (*evaluate_internal_energy)(
      PetscScalar* V_system,               //!
      const Vec mean_q,                    //!
      const Vec xi,                        //!
      const Vec rho,                       //!
      const AtomicSpecie* specie_ptr,      //!
      const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute \f$\frac{\partial V}{\partial \mathbf{q}_{i^*}}\f$:
   * \f[
   *  \frac{\partial V}{\partial
   * \mathbf{q}_{i^*}}\ =\ \frac{\partial}{\partial \mathbf{q}_{i^*}} \sum_i
   * V_i \f] Derivative of the total potential with respect the site
   * position
   * @param mean_q Position (\f$\mathbf{q}\f$) of the atomic
   * positions in the updated configuration
   * @param xi Molar fraction
   * @param rho Local density of the system
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return RHS Vector with the derivative of the potential
   */
  PetscErrorCode (*evaluate_dV_dq)(Vec RHS,                             //!
                                   const Vec mean_q,                    //!
                                   const Vec xi,                        //!
                                   const Vec rho,                       //!
                                   const AtomicSpecie* specie_ptr,      //!
                                   const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute \f$\frac{\partial^2 V}{\partial \mathbf{q}_{i^*}^2}\f$:
   * \f[
   *  \frac{\partial^2 V}{\partial
   * \mathbf{q}_{i^*}^2}\ =\ \frac{\partial}{\partial \mathbf{q}_{i^*}} \sum_i
   * \frac{\partial V_i}{\partial \mathbf{q}_{i^*}} \f]
   * Derivative of the total potential with respect the site position
   * @param mean_q Position (\f$\mathbf{q}\f$) of the atomic
   * positions in the updated configuration
   * @param xi Molar fraction
   * @param rho Local density of the system
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return Jac Jacobian matrix with the second derivative of the potential
   */
  PetscErrorCode (*evaluate_d2V_dq_ii)(Mat Jac,                             //!
                                       const Vec mean_q,                    //!
                                       const Vec xi,                        //!
                                       const Vec rho,                       //!
                                       const AtomicSpecie* specie_ptr,      //!
                                       const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the standard desviation of the position at each site using a
   * quasi-harmonic approximation of the thermalized potential
   * @param mean_q Position (\f$\mathbf{q}\f$) of the atomic
   * positions in the updated configuration
   * @param xi Molar fraction
   * @param rho Local density of the system
   * @param beta Lagrange multiplier (thermal)
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return stdv_q Standard deviation of the site position
   */
  PetscErrorCode (*evaluate_qh_stdv_q)(Vec stdv_q,                          //!
                                       const Vec mean_q,                    //!
                                       const Vec xi,                        //!
                                       const Vec rho,                       //!
                                       const Vec beta,                      //!
                                       const AtomicSpecie* specie_ptr,      //!
                                       const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the derivative of the potential functional with respect
   * the right stretch tensor:
   * \f[
   * \frac{\partial V}{\partial \mathbf{F}}
   * \f]
   * @param mean_q_ref Position (\f$\mathbf{q}_0\f$) of the atomic
   * positions in the reference configuration
   * @param mean_q Position (\f$\mathbf{q}\f$) of the atomic
   * positions in the updated configuration
   * @param xi Molar fraction
   * @param rho Local density of the system
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return RHS Vector with the derivative of the potential with respect to F
   */
  PetscErrorCode (*evaluate_dV_dF)(Vec RHS,                             //!
                                   const Vec mean_q_ref,                //!
                                   const Vec mean_q,                    //!
                                   const Vec xi,                        //!
                                   const Vec rho,                       //!
                                   const AtomicSpecie* specie_ptr,      //!
                                   const AtomTopology* atom_topology);  //!

  /**
   * @brief Evaluate the meanfield local density of the system
   * * @param mean_q:Mean value of q
   * @param stdv_q Standard desviation of q
   * @param xi Molar fraction
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return rho
   */
  PetscErrorCode (*evaluate_mf_rho)(Vec mf_rho,                          //!
                                    const Vec mean_q,                    //!
                                    const Vec stdv_q,                    //!
                                    const Vec xi,                        //!
                                    const AtomicSpecie* specie_ptr,      //!
                                    const AtomTopology* atom_topology);  //!

  /** @brief Total thermalized potential of the system
   * @param V0_system Total meanfield potential of the system
   * @param mean_q Mean value of q
   * @param stdv_q Standard desviation of q
   * @param xi Molar fraction
   * @param mf_rho Meanfield energy density
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return V0_system Total meanfield potential of the system
   */
  PetscErrorCode (*evaluate_mf_internal_energy)(
      PetscScalar* V0_system,              //!
      const Vec mean_q,                    //!
      const Vec stdv_q,                    //!
      const Vec xi,                        //!
      const Vec mf_rho,                    //!
      const AtomicSpecie* specie_ptr,      //!
      const AtomTopology* atom_topology);  //!

  /** @brief Evaluate the stress tensor of the system
   * @param Stress Stress tensor of the system
   * @param mean_q Mean value of q
   * @param stdv_q Standard desviation of q
   * @param xi Molar fraction
   * @param beta Lagrange multiplier (thermal)
   * @param mf_rho Meanfield energy density
   * @param specie_ptr Information of the atomic species
   * @param atom_topology Topology information of the atoms
   * @return Stress Stress tensor of the system
   */
  PetscErrorCode (*evaluate_virial_stress)(
      Vec Stress,                          //!
      const Vec mean_q,                    //!
      const Vec stdv_q,                    //!
      const Vec xi,                        //!
      const Vec mf_rho,                    //!
      const AtomicSpecie* specie_ptr,      //!
      const AtomTopology* atom_topology);  //!

  /**
   * @brief Evaluate the meanfield entropy of the system
   * @param S0_system Meanfield entropy of the system
   * @param mean_q Mean value of q
   * @param stdv_q Standard desviation of q
   * @param xi Molar fraction
   * @param beta Lagrange multiplier (thermal)
   * @param gamma Lagrange multiplier (chemical)
   * @param specie_ptr Information of the atomic species
   * @return S0_system Meanfield entropy of the system
   */
  PetscErrorCode (*evaluate_entropy)(PetscScalar* S0_system,           //!
                                     const Vec mean_q,                 //!
                                     const Vec stdv_q,                 //!
                                     const Vec xi,                     //!
                                     const Vec beta,                   //!
                                     const Vec gamma,                  //!
                                     const AtomicSpecie* specie_ptr);  //!

  /***
   * @brief Compute the derivatives of the thermalised functional with respect
   * the mean position of each atomis site i:
   *
   * \f[
   *    \frac{\partial V_0}{\partial \bar{\mathbf{q}}_i}
   * \f]
   *
   * @param RHS: Right-hand side of the derivative
   * @param mean_q: Mean value of q
   * @param stdv_q: Standard desviation of q
   * @param xi: Molar fraction
   * @param beta: Lagrange multiplier (thermal)
   * @param mf_rho: Meanfield energy density
   * @param specie_ptr: Information of the atomic species
   * @param atom_topology: Topology information of the atoms
   * @return RHS Vector with the derivative of the potential
   */
  PetscErrorCode (*evaluate_dL0_dmq)(Vec RHS,                             //!
                                     const Vec mean_q,                    //!
                                     const Vec stdv_q,                    //!
                                     const Vec xi,                        //!
                                     const Vec beta,                      //!
                                     const Vec mf_rho,                    //!
                                     const AtomicSpecie* specie_ptr,      //!
                                     const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the hessian of the potential with respect the mean position
   * of the site positions.
   *
   * \f[
   *    \frac{\partial^2 V_0}{\partial \bar{\mathbf{q}}_i^2}
   * \f]
   *
   * @param Jac: Jacobian matrix of the derivative
   * @param mean_q: Mean value of q
   * @param stdv_q: Standard desviation of q
   * @param xi: Molar fraction
   * @param mf_rho: Meanfield energy density
   * @param beta: Lagrange multiplier (thermal)
   * @param specie_ptr: Information of the atomic species
   * @param atom_topology: Topology information of the atoms
   * @return Jac Vector with the derivative of the potential
   */
  PetscErrorCode (*evaluate_d2L0_dmq_ii)(
      Mat Jac,                             //!
      const Vec mean_q,                    //!
      const Vec stdv_q,                    //!
      const Vec xi,                        //!
      const Vec mf_rho,                    //!
      const Vec beta,                      //!
      const AtomicSpecie* specie_ptr,      //!
      const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the derivatives of the free-entropy functional with respect
   the standard desviation of each atomis site i:
   *
   \f[
    \frac{\partial V_0}{\partial \sigma_i}
   \f]
   * @param RHS: Right-hand side of the derivative
   * @param mean_q: Mean value of q
   * @param stdv_q: Standard desviation of q
   * @param xi: Molar fraction
   * @param beta: Lagrange multiplier (thermal)
   * @param mf_rho: Meanfield energy density
   * @param specie_ptr: Information of the atomic species
   * @param atom_topology: Topology information of the atoms
   * @return RHS Vector with the derivative of the potential
   */
  PetscErrorCode (*evaluate_dL0_dsq)(Vec RHS,                             //!
                                     const Vec mean_q,                    //!
                                     const Vec stdv_q,                    //!
                                     const Vec xi,                        //!
                                     const Vec beta,                      //!
                                     const Vec mf_rho,                    //!
                                     const AtomicSpecie* specie_ptr,      //!
                                     const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the derivatives of the free-entropy functional with respect
   * the mean and standard desviation of the position at the atomic site i:
   *
   * \f[
   *    \frac{\partial V_0}{\partial \bar{\mathbf{q}}_i} \\
   *    \frac{\partial V_0}{\partial \sigma_i}
   * \f]
   * @param RHS_mean_q: Right-hand side of the derivative with respect the mean
   * position
   * @param RHS_stdv_q: Right-hand side of the derivative with respect the
   * standard desviation
   * @param mean_q: Mean value of q
   * @param stdv_q: Standard desviation of q
   * @param xi: Molar fraction
   * @param beta: Lagrange multiplier (thermal)
   * @param mf_rho: Meanfield energy density
   * @param specie_ptr: Information of the atomic species
   * @param atom_topology: Topology information of the atoms
   * @return RHS_mean_q Vector with the derivative of the potential with respect
   * to the mean position
   * @return RHS_stdv_q Vector with the derivative of the potential with respect
   * to the standard desviation
   */
  PetscErrorCode (*evaluate_dL0_dmq_dsq)(
      Vec RHS_mean_q,                      //!
      Vec RHS_stdv_q,                      //!
      const Vec mean_q,                    //!
      const Vec stdv_q,                    //!
      const Vec xi,                        //!
      const Vec beta,                      //!
      const Vec mf_rho,                    //!
      const AtomicSpecie* specie_ptr,      //!
      const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the derivatives of the free-entropy functional with respect
   * the deformation gradient:
   *
   \f[
    \frac{\partial V_0}{\partial \mathbf{F}}
   \f]
   * @param RHS: Right-hand side of the derivative
   * @param mean_q: Mean value of q
   * @param mean_q_ref: Mean value of q in the reference configuration
   * @param stdv_q: Standard deviation of q
   * @param xi: Molar fraction
   * @param beta: Lagrange multiplier (thermal)
   * @param mf_rho: Meanfield energy density
   * @param specie_ptr: Information of the atomic species
   * @param atom_topology: Topology information of the atoms
   * @return RHS Vector with the derivative of the potential with respect to the
   deformation gradient
   */
  PetscErrorCode (*evaluate_dL0_dF)(Vec RHS,
                                    const Vec mean_q,                    //!
                                    const Vec mean_q_ref,                //!
                                    const Vec stdv_q,                    //!
                                    const Vec xi,                        //!
                                    const Vec beta,                      //!
                                    const Vec mf_rho,                    //!
                                    const AtomicSpecie* specie_ptr,      //!
                                    const AtomTopology* atom_topology);  //!

  /**
   * @brief
   *
   */
  PetscErrorCode (*evaluate_free_entropy)(
      PetscScalar* free_S0_system,         //!
      const Vec mean_q,                    //!
      const Vec stdv_q,                    //!
      const Vec xi,                        //!
      const Vec mf_rho,                    //!
      const Vec beta,                      //!
      const Vec gamma,                     //!
      const AtomicSpecie* specie_ptr,      //!
      const AtomTopology* atom_topology);  //!

  PetscErrorCode (*evaluate_dxdt)(Vec dxdt,                           //!
                                  const Vec mean_q,                   //!
                                  const Vec stdv_q,                   //!
                                  const Vec xi,                       //!
                                  const Vec xi_ref,                   //!
                                  const Vec beta,                     //!
                                  const Vec gamma,                    //!
                                  const PetscInt* gamma_bcc_ptr,      //!
                                  const PetscInt* idx_diff_ptr,       //!
                                  const AtomTopology* atom_topology,  //!
                                  double* particle_mass);             //!

  PetscErrorCode (*evaluate_Jac_dxdt)(Mat Jac,                            //!
                                      const Vec mean_q,                   //!
                                      const Vec stdv_q,                   //!
                                      const Vec xi,                       //!
                                      const Vec xi_ref,                   //!
                                      const Vec beta,                     //!
                                      const Vec gamma,                    //!
                                      const PetscInt* gamma_bcc_ptr,      //!
                                      const PetscInt* idx_diff_ptr,       //!
                                      const AtomTopology* atom_topology,  //!
                                      double* particle_mass);             //!

  /**
   * @brief Compute the discrete Jordan-Kinderlerhrer-Otto functional
   *
   \f[
    F(\rho_{k+1}) \sim \sum_p^n \frac{1}{2} \frac{\|q_{p,k+1} - q_{p,k} \|^2 +
   \|\beta_{p,k+1} - \beta_{p,k} \|^2}{t_{k+1} - t_{k}} m_p +
   V(\rho_{k+1}(x_{p,k+1}))\, m_p
   \f]
   * @param JKO_system: Value of the JKO functional
   * @param Delta_t: Time step increment t_{k+1} - t_{k}
   * @param rho: Density field
   * @param x_k1: Particle position at k+1
   * @param x_k: Particle position at k
   * @param beta_k1: Particle width at k+1
   * @param beta_k: Particle width at k
   * @param atom_topology: Topology information of the atoms
   */
  PetscErrorCode (*evaluate_meassure_JKO)(
      Vec rho,                             //!
      const Vec q_k1,                      //!
      const Vec beta_k1,                   //!
      const Vec mass,                      //!
      const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the discrete Jordan-Kinderlerhrer-Otto functional
   *
   \f[
    F(\rho_{k+1}) \sim \sum_p^n \frac{1}{2} \frac{\|q_{p,k+1} - q_{p,k} \|^2 +
   \|\beta_{p,k+1} - \beta_{p,k} \|^2}{t_{k+1} - t_{k}} m_p +
   V(\rho_{k+1}(x_{p,k+1}))\, m_p
   \f]
   * @param JKO_system: Value of the JKO functional
   * @param rho: Density field
   * @param x_k1: Particle position at k+1
   * @param x_k: Particle position at k
   * @param beta_k1: Particle width at k+1
   * @param beta_k: Particle width at k
   * @param atom_topology: Topology information of the atoms
   */
  PetscErrorCode (*evaluate_JKO)(PetscScalar* JKO_system,             //!
                                 PetscScalar Delta_t,                 //!
                                 const Vec rho,                       //!
                                 const Vec q_k1,                      //!
                                 const Vec q_k,                       //!
                                 const Vec beta_k1,                   //!
                                 const Vec beta_k,                    //!
                                 const Vec mass,                      //!
                                 const AtomTopology* atom_topology);  //!

  /**
   * @brief Compute the discrete Jordan-Kinderlerhrer-Otto functional
   *
   \f[
    F(\rho_{k+1}) \sim \sum_p^n \frac{1}{2} \frac{\|q_{p,k+1} - q_{p,k} \|^2 +
   \|\beta_{p,k+1} - \beta_{p,k} \|^2}{t_{k+1} - t_{k}} m_p +
   V(\rho_{k+1}(q_{p,k+1}))\, m_p
   \f]
   * @param D_JKO_Dq: Value of the position gradient of the JKO functional
   * @param rho: Density field
   * @param x_k1: Particle position at k+1
   * @param x_k: Particle position at k
   * @param beta_k1: Particle width at k+1
   * @param atom_topology: Topology information of the atoms
   */
  PetscErrorCode (*evaluate_D_JKO_Dq)(
      Vec D_JKO_Dq,                        //!
      PetscScalar Delta_t,                 //!
      const Vec rho_k1,                    //!
      const Vec x_k1,                      //!
      const Vec x_k,                       //!
      const Vec beta_k1,                   //!
      const Vec mass,                      //!
      const AtomTopology* atom_topology);  //!                  //!

} dmd_equations;

/*******************************************************/

/**
 * @brief Abstract base class for Shape Functions and Potential Components.
 * * This interface defines the contract for any function N_i and its
 * spatial gradient DN_i, evaluated at point x relative to a reference
 * position x_i (e.g., a node or an atomic site).
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
  virtual double N_i(const Eigen::Vector3d& x, const Eigen::Vector3d& x_i,
                     double beta_i) const = 0;

  /**
   * @brief Calculates the gradient (first derivative) of the function.
   * * @param x The evaluation point coordinates.
   * @param x_i The reference position or nodal coordinates.
   * @return A Vector3d containing the partial derivatives [dN/dx, dN/dy,
   * dN/dz].
   */
  virtual Eigen::Vector3d dN_i(const Eigen::Vector3d& x,
                               const Eigen::Vector3d& x_i,
                               double beta_i) const = 0;
};

/*******************************************************/

typedef struct DMD_context {

  //! @param n_sites_global: Number of atoms in the global domain
  PetscInt n_sites_global;

  //! @param n_sites_local: Number of atoms in the local domain (without ghost)
  PetscInt n_sites_local;

  //! @param n_ghost: Number of ghost atoms in the local domain
  PetscInt n_ghost;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  @brief Interatomic Potential Information
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  PetscScalar r_cutoff;

  //  adpPotential V;

  dmd_equations system_equations;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  @brief Enviroment thermo-chemo-mechanical variables
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  //! @param Pressure_env: Enviromental value of the pressure
  double Pressure_env;

  //! @param Temperature_env: Enviromental value of the temperature
  double Temperature_env;

  //! @param ChemicalPotential_env: Enviromental value of the chemical potential
  //! ({mu}). {mu} = {gamma}/{beta}
  double ChemicalPotential_env;

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @brief Topology variables
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  /** @param dump2petsc_mapping: Allow the user to track the index of each atom
   * from the .dump prdering and petsc ordering */
  AO dump2petsc_mapping;

} DMD_context;

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
#define r_cutoff_ADP_MgHx 6.29867867868
#define r_cutoff_Angular_ADP_MgHx 4.0
#define r_cutoff_ADP_AlCu 6.2872100000
#define r_cutoff_EAM_FS_FeHx 5.3001060021200423  //
#define r_cutoff_Eb 3.0                          // 3.2 //

//! Virtual indenter (rigid body) repulsive pair parameters
#define r_cutoff_Indenter_default 4.5
#define sigma_Indenter_default 1.5
#define epsilon_Indenter_default 0.1

#define maxneigh 500
#define max_chemical_neighs 250

//! This values has been computed in a MgHx hcp cell and they
//! are the values which warraty the stability of the system. If the occupancy
//! is below this value do not compute the equilibrium Zero value for the
//! occupancy
#define max_occupancy 0.9999
#define min_occupancy 0.0001
#define min_occupancy_mean_q 1e-4
#define min_occupancy_stdv_q 1e-2
#define min_occupancy_diff 1e-3
#define min_stdv_q 0.001
#define max_stdv_q 3.0

// Define default min value for the energy barrier eV
#define min_Energy_barrier_default 0.2

/*
  Math macros
*/
static float sqr_arg;
#define SQR(a) ((sqr_arg = (a)) == 0.0 ? 0.0 : sqr_arg * sqr_arg)
static double dsqr_arg;
#define DSQR(a) ((dsqr_arg = (a)) == 0.0 ? 0.0 : dsqr_arg * dsqr_arg)
static double dmax_arg1, dmax_arg2;
#define DMAX(a, b)                   \
  (dmax_arg1 = (a), dmax_arg2 = (b), \
   (dmax_arg1) > (dmax_arg2) ? (dmax_arg1) : (dmax_arg2))
static double dmin_arg1, dmin_arg2;
#define DMIN(a, b)                   \
  (dmin_arg1 = (a), dmin_arg2 = (b), \
   (dmin_arg1) < (dmin_arg2) ? (dmin_arg1) : (dmin_arg2))
static float max_arg1, max_arg2;
#define FMAX(a, b)                 \
  (max_arg1 = (a), max_arg2 = (b), \
   (max_arg1) > (max_arg2) ? (max_arg1) : (max_arg2))
static float min_arg1, min_arg2;
#define FMIN(a, b)                 \
  (min_arg1 = (a), min_arg2 = (b), \
   (min_arg1) < (min_arg2) ? (min_arg1) : (min_arg2))
static long lmax_arg1, lmax_arg2;
#define LMAX(a, b)                   \
  (lmax_arg1 = (a), lmax_arg2 = (b), \
   (lmax_arg1) > (lmax_arg2) ? (lmax_arg1) : (lmax_arg2))
static long lmin_arg1, lmin_arg2;
#define LMIN(a, b)                   \
  (lmin_arg1 = (a), lmin_arg2 = (b), \
   (lmin_arg1) < (lmin_arg2) ? (lmin_arg1) : (lmin_arg2))
static int imax_arg1, imax_arg2;
#define IMAX(a, b)                   \
  (imax_arg1 = (a), imax_arg2 = (b), \
   (imax_arg1) > (imax_arg2) ? (imax_arg1) : (imax_arg2))
static int imin_arg1, imin_arg2;
#define IMIN(a, b)                   \
  (imin_arg1 = (a), imin_arg2 = (b), \
   (imin_arg1) < (imin_arg2) ? (imin_arg1) : (imin_arg2))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) s

#endif
/**
 * @file Potentials/Advection-Diff-OpenMP.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief Advection–diffusion realization of GoverningEquations.
 * @version 0.1
 * @date 2023-05-19
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef ADVEC_DIFF_POTENTIAL_HPP
#define ADVEC_DIFF_POTENTIAL_HPP

#include "Macros.hpp"

/**
 * @brief Advection–diffusion JKO equations (entropy + Wasserstein transport).
 *
 * Concrete GoverningEquations used by the blob mass-transport solver. Derive
 * further classes from GoverningEquations (or from this one) to swap potentials
 * or transport models without changing the TAO driver.
 */
class AdvectionDiffusionEquations : public GoverningEquations {
 public:
  ~AdvectionDiffusionEquations() override {}

  PetscErrorCode evaluate_meassure_JKO(
      Vec rho, const Vec q_k1, const Vec beta_k1, const Vec mass,
      const ParticleTopology* atom_topology) override;

  PetscErrorCode evaluate_JKO(PetscScalar* JKO_system, PetscScalar Delta_t,
                              const Vec rho, const Vec q_k1, const Vec q_k,
                              const Vec beta_k1, const Vec beta_k,
                              const Vec mass,
                              const ParticleTopology* atom_topology) override;

  PetscErrorCode evaluate_D_JKO_Dq(
      Vec D_JKO_Dq, PetscScalar Delta_t, const Vec rho_k1, const Vec x_k1,
      const Vec x_k, const Vec beta_k1, const Vec mass,
      const ParticleTopology* atom_topology) override;
};

#endif /* ADVEC_DIFF_POTENTIAL_HPP */

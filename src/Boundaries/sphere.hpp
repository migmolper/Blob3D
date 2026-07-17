/**
 * @file Boundaries/sphere.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief Sphere boundary condition.
 * @version 0.1
 * @date 2023-05-19
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef BC_SPHERE_HPP
#define BC_SPHERE_HPP

#include "Macros.hpp"

class BC_Sphere : public boundaryCondition {
public:
  /**
   * @brief Constructor of the Sphere class
   * @param radius_sphere Radius of the sphere
   * @param center_sphere Center of the sphere
   * @param penalty Penalty stiffness
   * @param buffer Buffer distance
   */
  BC_Sphere(PetscScalar radius_sphere, const Eigen::Vector3d center_sphere,
            PetscScalar penalty, PetscScalar buffer)
      : radius_sphere(radius_sphere), center_sphere(center_sphere),
        penalty(penalty), buffer(buffer) {}

  /**
   * @brief Destructor
   */
  ~BC_Sphere() override {}

  /**
   * @brief Add the barrier potential to the JKO system
   * @param JKO_system JKO system
   * @param x_k1 Position of the particles at the previous time step
   * @param mass Mass of the particles
   * @return PetscErrorCode
   */
  PetscErrorCode
  add_barrier_potential(PetscScalar *JKO_system, const Vec x_k1, const Vec mass) override;

  /**
   * @brief Add the barrier forces to the JKO system
   * @param D_JKO_Dq Derivative of the JKO system
   * @param x_k1 Position of the particles at the previous time step
   * @param mass Mass of the particles
   * @return PetscErrorCode
   */
  PetscErrorCode
  add_barrier_forces(Vec D_JKO_Dq, const Vec x_k1, const Vec mass) override;

public:
  /** @brief Radius of the sphere */
  PetscScalar radius_sphere;
  /** @brief Penalty stiffness */
  PetscScalar penalty;
  /** @brief Buffer distance */
  PetscScalar buffer;
  /** @brief Center of the sphere */
  const Eigen::Vector3d center_sphere;
};

#endif /* SPHERE_HPP */
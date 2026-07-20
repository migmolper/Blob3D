/**
 * @file sphere.hpp
 * @author Miguel Molinos (@migmolper)
 * @brief Sphere boundary condition.
 * @version 0.1
 * @date 2023-05-19
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "Macros.hpp"

class BC_Sphere : public boundaryCondition {
public:
  /**
   * @brief Constructor of the Sphere class
   * @param radius_sphere Radius of the sphere
   * @param center_sphere Center of the sphere
   * @param penalty Penalty stiffness
   */
  BC_Sphere(PetscScalar radius_sphere, const Eigen::Vector3d center_sphere,
            PetscScalar penalty)
      : radius_sphere(radius_sphere), center_sphere(center_sphere),
        penalty(penalty) {}

  /**
   * @brief Destructor
   */
  ~BC_Sphere() override {}

  /**
   * @brief Calculate the potential between the particle and the boundary.
   * @param x The particle position.
   * @return The potential between the particle and the boundary.
   */
  PetscScalar potential(const Eigen::Vector3d &x) const override {
    return 0.5 * penalty * DSQR(this->distance(x));
  }

  /**
   * @brief Calculate the gradient of the potential between the particle and the
   * boundary.
   *
   * @param x The particle position.
   * @return The gradient of the potential between the particle and the
   * boundary.
   */
  Eigen::Vector3d gradient_potential(const Eigen::Vector3d &x) const override {
    return penalty * this->distance(x) * this->gradient_distance(x);
  }

  /** @brief Radius of the sphere */
  PetscScalar radius_sphere;
  /** @brief Penalty stiffness */
  PetscScalar penalty;
  /** @brief Center of the sphere */
  const Eigen::Vector3d center_sphere;

private:
  /**
   * @brief Calculate the distance between the particle and the boundary.
   * @param x The particle position.
   * @return The distance between the particle and the boundary.
   */
  double distance(const Eigen::Vector3d &x) const {
    return radius_sphere - (x - center_sphere).norm();
  }

  /**
   * @brief Calculate the gradient of the distance between the particle and the
   * boundary.
   *
   * @param x The particle position.
   * @return The gradient of the distance between the particle and the boundary.
   */
  Eigen::Vector3d gradient_distance(const Eigen::Vector3d &x) const {
    return - (x - center_sphere).normalized();
  }
};

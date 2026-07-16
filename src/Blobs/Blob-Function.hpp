/**
 * @file Blob-Function.hpp
 * @brief Gaussian Blob implementation derived from ShapeFunction.
 * @date 2026-03-21
 */

#ifndef BLOB_FUNCTION_HPP
#define BLOB_FUNCTION_HPP

#include <Eigen/Dense>
#include "Macros.hpp" 

class Blob : public ShapeFunction {
public:
    virtual ~Blob() {}

    /**
     * @brief Implements N_i(x, x_i) = norm * exp(-beta * |x - x_i|^2)
     */
    double N_i(const Eigen::Vector3d& x, const Eigen::Vector3d& x_i,
               double beta_i) const override;

    /**
     * @brief Implements DN_i = -2 * beta * N_i * (x - x_i)
     */
    Eigen::Vector3d dN_i(const Eigen::Vector3d& x, const Eigen::Vector3d& x_i,
                         double beta_i) const override;
};

#endif // BLOB_FUNCTION_HPP
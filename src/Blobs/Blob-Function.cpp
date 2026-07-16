/**
 * @file Blob-Function.cpp
 * @date 2026-03-21
 */

#include "Blobs/Blob-Function.hpp"
#include <cmath>

/************************************************************************/

double Blob::N_i(const Eigen::Vector3d& x, const Eigen::Vector3d& x_i,
                 double beta_i) const {

  double r2 = (x - x_i).squaredNorm();

  return std::pow((beta_i / M_PI), 1.5) * std::exp(-beta_i * r2);
}

/************************************************************************/

Eigen::Vector3d Blob::dN_i(const Eigen::Vector3d& x, const Eigen::Vector3d& x_i,
                           double beta_i) const {

  Eigen::Vector3d dx = x - x_i;

  return -N_i(x, x_i, beta_i) * 2.0 * beta_i * dx;
}

/************************************************************************/
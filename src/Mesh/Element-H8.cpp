/**
 * @file Element-H8.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2026-05-10
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "Macros.hpp"
#include <Eigen/Dense>
#include <petscdmplex.h>

/*******************************************************/

bool InOut_Element_H8(const PetscScalar* mean_q_i,
                      const PetscScalar* el_coords) {

  PetscReal ll_vertex[3], ur_vertex[3];

  ll_vertex[0] = el_coords[0];
  ll_vertex[1] = el_coords[1];
  ll_vertex[2] = el_coords[2];
  ur_vertex[0] = el_coords[0];
  ur_vertex[1] = el_coords[1];
  ur_vertex[2] = el_coords[2];

  for (int i = 1; i < 8; i++) {
    ll_vertex[0] = DMIN(ll_vertex[0], el_coords[i * 3 + 0]);
    ll_vertex[1] = DMIN(ll_vertex[1], el_coords[i * 3 + 1]);
    ll_vertex[2] = DMIN(ll_vertex[2], el_coords[i * 3 + 2]);
    ur_vertex[0] = DMAX(ur_vertex[0], el_coords[i * 3 + 0]);
    ur_vertex[1] = DMAX(ur_vertex[1], el_coords[i * 3 + 1]);
    ur_vertex[2] = DMAX(ur_vertex[2], el_coords[i * 3 + 2]);
  }

  if ((mean_q_i[0] >= ll_vertex[0]) &&  //!
      (mean_q_i[0] < ur_vertex[0]) &&   //!
      (mean_q_i[1] >= ll_vertex[1]) &&  //!
      (mean_q_i[1] < ur_vertex[1]) &&   //!
      (mean_q_i[2] >= ll_vertex[2]) &&  //!
      (mean_q_i[2] < ur_vertex[2])) {
    return PETSC_TRUE;
  } else {
    return PETSC_FALSE;
  }
}

/*******************************************************/
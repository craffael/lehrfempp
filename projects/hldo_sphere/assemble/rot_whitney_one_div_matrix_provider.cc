#include "rot_whitney_one_div_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd RotWhitneyOneDivMatrixProvider::Eval(
    const lf::mesh::Entity &entity) const {
  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // correct for orientation
  auto edgeOrientations = entity.RelativeOrientations();
  double s0 = lf::mesh::to_sign(edgeOrientations[0]);
  double s1 = lf::mesh::to_sign(edgeOrientations[1]);
  double s2 = lf::mesh::to_sign(edgeOrientations[2]);

  // Compute the element matrix for the product of baricentric functions
  Eigen::MatrixXd elem_mat(3, 1);
  elem_mat << -s0, -s1, -s2;

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble

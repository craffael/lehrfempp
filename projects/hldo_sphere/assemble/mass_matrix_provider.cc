#include "mass_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd MassMatrixProvider::Eval(const lf::mesh::Entity &entity) {
  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get the geometry of the entity
  const auto *geom = entity.Geometry();

  // Compute the element matrix for the product of baricentric functions
  Eigen::MatrixXd elem_mat(3, 3);
  // clang-format off
  elem_mat <<   2, 1, 1,
                1, 2, 1,
                1, 1, 2;
  // clang-format on
  elem_mat *= lf::geometry::Volume(*geom) / 12.0;

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble

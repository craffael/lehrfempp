#include "whitney_two_mass_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd WhitneyTwoMassMatrixProvider::Eval(
    const lf::mesh::Entity &entity) const {
  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get the geometry of the entity
  const auto *geom = entity.Geometry();

  // Compute the element matrix for the product of baricentric functions
  Eigen::MatrixXd elem_mat(1, 1);
  // clang-format off
  elem_mat <<   1;
  // clang-format on
  elem_mat *= lf::geometry::Volume(*geom);

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble

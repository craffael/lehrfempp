#include "whitney_one_form_curl_element_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd WhitneyOneFormCurlElementMatrixProvider::Eval(
    const lf::mesh::Entity &entity) const {
  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get the geometry of the entity
  const auto *geom = entity.Geometry();

  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

  // fill the element matrix
  Eigen::MatrixXd elem_mat = Eigen::MatrixXd::Ones(3, 3);
  elem_mat *= 1 / lf::geometry::Volume(*geom);

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble

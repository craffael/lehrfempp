#include "whitney_one_curl_curl_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd WhitneyOneCurlCurlMatrixProvider::Eval(
    const lf::mesh::Entity &entity) const {
  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get the geometry of the entity
  const auto *geom = entity.Geometry();

  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

  // compute the orientations
  auto edgeOrientations = entity.RelativeOrientations();
  Eigen::VectorXd s(3);
  for (int i = 0; i < 3; i++) {
    s(i) = lf::mesh::to_sign(edgeOrientations[i]);
  }

  // fill the element matrix
  Eigen::MatrixXd elem_mat = Eigen::MatrixXd::Ones(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      elem_mat(i, j) = s(i) * s(j);
    }
  }

  elem_mat *= 1 / lf::geometry::Volume(*geom);

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble

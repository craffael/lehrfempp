#include "utils.h"

#include <lf/geometry/geometry.h>

namespace projects::ipdg_stokes::mesh {

Eigen::Matrix<double, 2, 3> computeOutwardNormals(
    const lf::mesh::Entity &entity) {
  // Get the geometry of the entity
  const auto geom = entity.Geometry();
  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());
  // Use the vertex coordinates to compute the local normals on the edges
  Eigen::Matrix<double, 2, 3> normals(vertices.rows(), vertices.cols());
  for (long vert = 0; vert < vertices.cols() - 1; ++vert) {
    normals(0, vert) = vertices(1, vert + 1) - vertices(1, vert + 0);
    normals(1, vert) = vertices(0, vert + 0) - vertices(0, vert + 1);
  }
  normals(0, vertices.cols() - 1) =
      vertices(1, 0) - vertices(1, vertices.cols() - 1);
  normals(1, vertices.cols() - 1) =
      vertices(0, vertices.cols() - 1) - vertices(0, 0);
  // Compute test vectors to test whether the normal faces inward or outward
  Eigen::Matrix<double, 2, 3> test;
  test << vertices.col(2) - vertices.col(0), vertices.col(0) - vertices.col(1),
      vertices.col(1) - vertices.col(2);
  // Compute the sign to flip the normals
  Eigen::Matrix<double, 1, 3> flip =
      (normals.array() * test.array()).matrix().colwise().sum();
  normals.array().row(0) *= flip.array();
  normals.array().row(1) *= flip.array();
  normals = -normals.colwise().normalized();
  return normals;
}

Eigen::Matrix<double, 2, 3> computeTangentials(const lf::mesh::Entity &entity) {
  // Compute the tangentials from the rotated normals
  const auto normals = computeOutwardNormals(entity);
  Eigen::Matrix<double, 2, 3> tang;
  tang << -normals.row(1), normals.row(0);
  return tang;
}

}  // end namespace projects::ipdg_stokes::mesh

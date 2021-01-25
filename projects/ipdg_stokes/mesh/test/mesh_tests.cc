#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <lf/base/span.h>
#include <lf/geometry/geometry.h>
#include <lf/geometry/tria_o1.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh.h>
#include <utils.h>

#include <array>

std::shared_ptr<lf::mesh::Mesh> buildTriangle(
    const Eigen::Matrix<double, 2, 3> &vertices) {
  lf::mesh::hybrid2d::MeshFactory factory(2);
  lf::base::size_type v_idx[3];
  v_idx[0] = factory.AddPoint(vertices.col(0));
  v_idx[1] = factory.AddPoint(vertices.col(1));
  v_idx[2] = factory.AddPoint(vertices.col(2));
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  factory.AddEntity(lf::base::RefEl::kTria(), v_idx, std::move(geom));
  return factory.Build();
}

TEST(projects_ipdg_stokes_mesh, triangle_counterclockwise) {
  // Obtain a reference triangle with the vertices numbered counterclockwise
  // starting at (0, 0)
  Eigen::Matrix<double, 2, 3> vertices;
  vertices << 0, 1, 0, 0, 0, 1;
  const auto mesh = buildTriangle(vertices);
  const auto entity = mesh->EntityByIndex(0, 0);
  // Compute the outward pointing normals
  const auto normals =
      projects::ipdg_stokes::mesh::computeOutwardNormals(*entity);
  // These are the reference normals
  Eigen::Matrix<double, 2, 3> normals_ref;
  normals_ref << 0, 1. / std::sqrt(2), -1, -1, 1. / std::sqrt(2), 0;
  // Test for equivalence
  EXPECT_EQ(normals.rows(), normals_ref.rows());
  EXPECT_EQ(normals.cols(), normals_ref.cols());
  for (int i = 0; i < normals_ref.rows(); ++i)
    for (int j = 0; j < normals_ref.cols(); ++j)
      EXPECT_DOUBLE_EQ(normals(i, j), normals_ref(i, j))
          << "mismatch in normals of edge " << j << " : ["
          << normals.col(j).transpose() << "] != ["
          << normals_ref.col(j).transpose() << "]";
}

TEST(projects_ipdg_stokes_mesh, triangle_clockwise) {
  // Obtain a reference triangle with the vertices numbered counterclockwise
  // starting at (0, 0)
  Eigen::Matrix<double, 2, 3> vertices;
  vertices << 0, 0, 1, 0, 1, 0;
  const auto mesh = buildTriangle(vertices);
  const auto entity = mesh->EntityByIndex(0, 0);
  // Compute the outward pointing normals
  const auto normals =
      projects::ipdg_stokes::mesh::computeOutwardNormals(*entity);
  // These are the reference normals
  Eigen::Matrix<double, 2, 3> normals_ref;
  normals_ref << -1, 1. / std::sqrt(2), 0, 0, 1. / std::sqrt(2), -1;
  // Test for equivalence
  EXPECT_EQ(normals.rows(), normals_ref.rows());
  EXPECT_EQ(normals.cols(), normals_ref.cols());
  for (int i = 0; i < normals_ref.rows(); ++i)
    for (int j = 0; j < normals_ref.cols(); ++j)
      EXPECT_DOUBLE_EQ(normals(i, j), normals_ref(i, j))
          << "mismatch in normals of edge " << j << " : ["
          << normals.col(j).transpose() << "] != ["
          << normals_ref.col(j).transpose() << "]";
}

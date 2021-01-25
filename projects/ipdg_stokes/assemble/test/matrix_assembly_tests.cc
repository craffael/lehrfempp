#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/base/base.h>
#include <lf/base/span.h>
#include <lf/geometry/geometry.h>
#include <lf/geometry/tria_o1.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/all_codim_mesh_data_set.h>
#include <lf/mesh/utils/utils.h>
#include <piecewise_const_element_matrix_provider.h>

#include <array>
#include <cmath>

TEST(projects_ipdg_stokes_assembly, piecewise_const_matrix_assembler_test) {
  // Build a mesh containing only the reference triangle
  const auto trig = lf::base::RefEl::kTria();
  const auto vertices = trig.NodeCoords();
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  lf::mesh::hybrid2d::MeshFactory factory(2);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const lf::mesh::utils::AllCodimMeshDataSet<bool> boundary(
      mesh, false);  // Assume that this triangle does not touch the boundary
  const double sigma = 1;
  const auto elem_mat_provider =
      projects::ipdg_stokes::assemble::PiecewiseConstElementMatrixProvider(
          sigma, boundary);
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(6, 6);
  const double r = std::sqrt(2);
  const double b = 1. / std::sqrt(2);
  const double s = 1. / (2 * sigma);
  // clang-format off
  Ae_anal <<  0,  0,  0, -1,  r, -1,
              0,  0,  0,  0, -b,  1,
	      0,  0,  0,  1, -b,  0,
	     -1,  0,  1, -s,  0,  0,
	      r, -b, -b,  0, -s,  0,
	     -1,  1,  0,  0,  0, -s;
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(Ae(i, j), Ae_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

TEST(projects_ipdg_stokes_assembly, global_matrix_assembly_1) {
  // Building a mesh out of two triangles such that
  // the shared edge has different orientations for each of them
  lf::mesh::hybrid2d::MeshFactory factory(2);
  Eigen::Matrix<double, 2, 4> vertices;
  vertices << -1, 0, 1, 0, 0, -1, 0, 1;
  lf::base::size_type v_idx[4];
  // Add the vertices
  for (int i = 0; i < vertices.cols(); ++i)
    v_idx[i] = factory.AddPoint(vertices.col(i));
  // Add the edges
  auto add_edge = [&](unsigned v1, unsigned v2) {
    const std::array<lf::base::size_type, 2> indices = {v1, v2};
    Eigen::Matrix<double, 2, 2> verts;
    verts << vertices.col(v1), vertices.col(v2);
    std::unique_ptr<lf::geometry::Geometry> geom =
        std::make_unique<lf::geometry::SegmentO1>(verts);
    factory.AddEntity(lf::base::RefEl::kSegment(),
                      nonstd::span(indices.data(), 2), std::move(geom));
  };
  add_edge(0, 1);
  add_edge(1, 3);
  add_edge(3, 0);
  add_edge(1, 2);
  add_edge(2, 3);
  // Add the elements
  auto add_triangle = [&](unsigned v1, unsigned v2, unsigned v3) {
    const std::array<lf::base::size_type, 3> indices = {v1, v2, v3};
    Eigen::Matrix<double, 2, 3> verts;
    verts << vertices.col(v1), vertices.col(v2), vertices.col(v3);
    std::unique_ptr<lf::geometry::Geometry> geom =
        std::make_unique<lf::geometry::TriaO1>(verts);
    factory.AddEntity(lf::base::RefEl::kTria(), nonstd::span(indices.data(), 3),
                      std::move(geom));
  };
  add_triangle(0, 1, 3);
  add_triangle(2, 3, 1);
  const auto mesh = factory.Build();

  // Compute the system matrix for the mesh
  lf::assemble::UniformFEDofHandler dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1}, {lf::base::RefEl::kSegment(), 1}});
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  const double sigma = 10;
  projects::ipdg_stokes::assemble::PiecewiseConstElementMatrixProvider
      elem_mat_provider(sigma, boundary);
  lf::assemble::COOMatrix<double> A_coo(dofh.NumDofs(), dofh.NumDofs());
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elem_mat_provider, A_coo);
  const Eigen::MatrixXd A = A_coo.makeDense().unaryExpr(
      [](double x) { return (std::fabs(x) < 1e-100 ? 0. : x); });

  // Construct the analytical system matrix
  Eigen::MatrixXd A_anal(9, 9);
  const double r = 1. / std::sqrt(2);
  const double s = 1. / sigma;
  // clang-format off
  A_anal <<  0,  0,  0,  0, -r,  1, -r,  0,  0,
	     0,  0,  0,  0,  0, -1,  r,  0,  r,
	     0,  0,  0,  0,  0,  1,  0, -r, -r,
	     0,  0,  0,  0,  r, -1,  0,  r,  0,
	    -r,  0,  0,  r, -s,  0,  0,  0,  0,
	     1, -1,  1, -1,  0, -s,  0,  0,  0,
	    -r,  r,  0,  0,  0,  0, -s,  0,  0,
	     0,  0, -r,  r,  0,  0,  0, -s,  0,
	     0,  r, -r,  0,  0,  0,  0,  0, -s;
  // clang-format on
  // Assert that the two matrices are approximately equal
  std::cout << A << std::endl;
  ASSERT_EQ(A.rows(), A_anal.rows());
  ASSERT_EQ(A.cols(), A_anal.cols());
  for (int i = 0; i < A_anal.rows(); ++i) {
    for (int j = 0; j < A_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(A(i, j), A_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

TEST(projects_ipdg_stokes_assembly, global_matrix_assembly_2) {
  // Building a mesh out of two triangles such that
  // the shared edge has the same orientation for each of them
  lf::mesh::hybrid2d::MeshFactory factory(2);
  Eigen::Matrix<double, 2, 4> vertices;
  vertices << -1, 0, 1, 0, 0, -1, 0, 1;
  lf::base::size_type v_idx[4];
  // Add the vertices
  for (int i = 0; i < vertices.cols(); ++i)
    v_idx[i] = factory.AddPoint(vertices.col(i));
  // Add the edges
  auto add_edge = [&](unsigned v1, unsigned v2) {
    const std::array<lf::base::size_type, 2> indices = {v1, v2};
    Eigen::Matrix<double, 2, 2> verts;
    verts << vertices.col(v1), vertices.col(v2);
    std::unique_ptr<lf::geometry::Geometry> geom =
        std::make_unique<lf::geometry::SegmentO1>(verts);
    factory.AddEntity(lf::base::RefEl::kSegment(),
                      nonstd::span(indices.data(), 2), std::move(geom));
  };
  add_edge(0, 1);
  add_edge(1, 3);
  add_edge(3, 0);
  add_edge(1, 2);
  add_edge(2, 3);
  // Add the elements
  auto add_triangle = [&](unsigned v1, unsigned v2, unsigned v3) {
    const std::array<lf::base::size_type, 3> indices = {v1, v2, v3};
    Eigen::Matrix<double, 2, 3> verts;
    verts << vertices.col(v1), vertices.col(v2), vertices.col(v3);
    std::unique_ptr<lf::geometry::Geometry> geom =
        std::make_unique<lf::geometry::TriaO1>(verts);
    factory.AddEntity(lf::base::RefEl::kTria(), nonstd::span(indices.data(), 3),
                      std::move(geom));
  };
  add_triangle(0, 1, 3);
  add_triangle(2, 1, 3);
  const auto mesh = factory.Build();

  // Compute the system matrix for the mesh
  lf::assemble::UniformFEDofHandler dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1}, {lf::base::RefEl::kSegment(), 1}});
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  const double sigma = 10;
  projects::ipdg_stokes::assemble::PiecewiseConstElementMatrixProvider
      elem_mat_provider(sigma, boundary);
  lf::assemble::COOMatrix<double> A_coo(dofh.NumDofs(), dofh.NumDofs());
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elem_mat_provider, A_coo);
  const Eigen::MatrixXd A = A_coo.makeDense().unaryExpr(
      [](double x) { return (std::fabs(x) <= 1e-100 ? 0. : x); });

  // Construct the analytical system matrix
  Eigen::MatrixXd A_anal(9, 9);
  const double r = 1. / std::sqrt(2);
  const double s = 1. / sigma;
  const double p = 0.5 / sqrt(2);
  // clang-format off
  A_anal <<  0,  0,  0,  0, -r,  1, -r,  0,  0,
	     0,  0,  0,  0,  0, -1,  r,  0,  r,
	     0,  0,  0,  0,  0,  1,  0, -r, -r,
	     0,  0,  0,  0,  r, -1,  0,  r,  0,
	    -r,  0,  0,  r, -s,  0,  0,  0,  0,
	     1, -1,  1, -1,  0, -s,  0,  0,  0,
	    -r,  r,  0,  0,  0,  0, -s,  0,  0,
	     0,  0, -r,  r,  0,  0,  0, -s,  0,
	     0,  r, -r,  0,  0,  0,  0,  0, -s;
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(A.rows(), A_anal.rows());
  ASSERT_EQ(A.cols(), A_anal.cols());
  std::cout << A_anal << std::endl;
  for (int i = 0; i < A_anal.rows(); ++i) {
    for (int j = 0; j < A_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(A(i, j), A_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

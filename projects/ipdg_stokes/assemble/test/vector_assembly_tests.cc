#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/geometry/tria_o1.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/all_codim_mesh_data_set.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <piecewise_const_element_vector_provider.h>

#include <array>
#include <span>

TEST(projects_ipdg_stokes_assembly, piecewise_const_vector_assembler_test) {
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
  factory.AddEntity(trig, std::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);

  // Apply constant volumetric forces in the x direction
  auto f = [](const Eigen::Vector2d &x) -> Eigen::Vector2d {
    Eigen::Vector2d force;
    force << 1, 0;
    return force;
  };
  // The dirichlet data will be unit vectors tangential to the edge
  auto dirichlet = [](const lf::mesh::Entity &edge) -> Eigen::Vector2d {
    const auto geom = edge.Geometry();
    const Eigen::Matrix2d verts = geom->Global(edge.RefEl().NodeCoords());
    return (verts.col(1) - verts.col(0)).normalized();
  };

  // Use the mispoint quadrature rule
  lf::quad::QuadRule qr = lf::quad::make_TriaQR_MidpointRule();

  // Test for no boundary edges
  {
    lf::mesh::utils::AllCodimMeshDataSet<bool> boundary(mesh, false);
    lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> dirichlet_mds(mesh, 1);
    for (const auto *ep : mesh->Entities(1)) {
      dirichlet_mds(*ep) = dirichlet(*ep);
    }
    const auto elem_vec_provider =
        projects::ipdg_stokes::assemble::PiecewiseConstElementVectorProvider(
            1, f, qr, boundary, dirichlet_mds);
    Eigen::VectorXd rhs = elem_vec_provider.Eval(*element);
    Eigen::VectorXd rhs_anal(6);
    rhs_anal << -0.5, 0, 0.5, 0, 0, 0;
    // Assert that the two vectors are approximately equal
    ASSERT_EQ(rhs.size(), rhs_anal.size());
    for (int i = 0; i < rhs_anal.size(); ++i)
      EXPECT_DOUBLE_EQ(rhs[i], rhs_anal[i]) << "mismatch in entry " << i;
  }
}

#include <gtest/gtest.h>

#include <norms.h>

#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>

TEST(projects_ipdg_stokes_post_processing, l2_constant) {
  // Build a mesh on $[0,1] \times [0,1]$
  std::shared_ptr<lf::mesh::MeshFactory> factory =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(factory);
  builder.setBottomLeftCorner(0, 0);
  builder.setTopRightCorner(1, 1);
  builder.setNoXCells(10);
  builder.setNoYCells(10);
  auto mesh = builder.Build();
  // Compute the L2 norm of the function constant zero on the mesh
  const double L2_zero = projects::ipdg_stokes::post_processing::L2norm(
      mesh, [](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return Eigen::Vector2d::Zero();
      });
  const double L2_one = projects::ipdg_stokes::post_processing::L2norm(
      mesh, [](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return Eigen::Vector2d::Ones();
      });
  ASSERT_NEAR(L2_zero, 0, 1e-10);
  ASSERT_NEAR(L2_one, std::sqrt(2), 1e-10);
}

TEST(projects_ipdg_stokes_post_processing, dg_constant) {
  // Build a mesh on $[0,1] \times [0,1]$
  std::shared_ptr<lf::mesh::MeshFactory> factory =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(factory);
  builder.setBottomLeftCorner(0, 0);
  builder.setTopRightCorner(1, 1);
  builder.setNoXCells(10);
  builder.setNoYCells(10);
  auto mesh = builder.Build();
  // Compute the DG norm of the function constant zero on the mesh
  const double DG_zero = projects::ipdg_stokes::post_processing::DGnorm(
      mesh,
      [](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return Eigen::Vector2d::Zero();
      },
      [](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return Eigen::Matrix2d::Zero();
      });
  const double DG_one = projects::ipdg_stokes::post_processing::DGnorm(
      mesh,
      [](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return Eigen::Vector2d::Ones();
      },
      [](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return Eigen::Matrix2d::Zero();
      });
  ASSERT_NEAR(DG_zero, 0, 1e-10);
  ASSERT_NEAR(DG_one, 0, 1e-10);
}

TEST(projects_ipdg_stokes_post_processing, dg_jump) {
  // Build a mesh on $[0,1] \times [0,1]$
  const size_t nx = 10;
  const size_t ny = 10;
  std::shared_ptr<lf::mesh::MeshFactory> factory =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(factory);
  builder.setBottomLeftCorner(0, 0);
  builder.setTopRightCorner(1, 1);
  builder.setNoXCells(nx);
  builder.setNoYCells(ny);
  auto mesh = builder.Build();
  // The function is zero for $y \in [0, 0.5]$ and one otherwise
  const double norm = projects::ipdg_stokes::post_processing::DGnorm(
      mesh,
      [&](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return mesh->Index(entity) < nx * ny ? Eigen::Vector2d::Ones()
                                             : Eigen::Vector2d::Zero();
      },
      [](const lf::mesh::Entity &entity, const Eigen::Vector2d &x) {
        return Eigen::Matrix2d::Zero();
      });
  ASSERT_NEAR(norm, std::sqrt(2 * nx), 1e-10);
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

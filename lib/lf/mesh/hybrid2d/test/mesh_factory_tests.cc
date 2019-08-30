/**
 * @file
 * @brief Tests for the lf::mesh::hybrid2d::MeshFactory
 * @author Raffael Casagrande
 * @date   2018-07-01 01:15:55
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <Eigen/Eigen>
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "mesh_factory_test.h"

namespace lf::mesh::hybrid2d::test {

bool mesh_sanity_check(const lf::mesh::Mesh& mesh) {
  std::cout << "Mesh sanity: Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(mesh);

  std::cout << "Mesh sanity: Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(mesh);

  std::cout << "Mesh sanity: Checking geometry compatibility: " << std::flush;
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;
  auto fails = lf::mesh::test_utils::isWatertightMesh(mesh, false);
  EXPECT_EQ(fails.size(), 0) << "Inconsistent geometry!";
  if (fails.size() == 0) {
    std::cout << "consistent!" << std::endl;
  } else {
    std::cout << "INCONSISTENT!" << std::endl;
    for (auto& geo_errs : fails) {
      std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                << std::endl;
    }
    return false;
  }

  // Compute volumes
  double total_area = 0.0;
  for (const mesh::Entity& cell : mesh.Entities(0)) {
    const double vol = Volume(*cell.Geometry());
    std::cout << cell.RefEl().ToString() << ' ' << mesh.Index(cell)
              << ": volume = " << vol << std::endl;
    total_area += vol;
  }
  std::cout << ">>> Total area = " << total_area << std::endl;

  std::cout << "Mesh sanity: Writing MATLAB file" << std::endl;
  io::writeMatlab(mesh, "test_mesh.m");

  // Printing mesh information
  lf::geometry::Geometry::output_ctrl_ = 20;
  lf::mesh::utils::printinfo_ctrl = 100;
  utils::PrintInfo(mesh, std::cout);
  return true;
}

// Test for generating a mesh with reconstruction of edge information
TEST(lf_edge_create, MeshFactory_p) {
  using coord_t = Eigen::Vector2d;
  using size_type = mesh::Mesh::size_type;

  // Building the mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  EXPECT_TRUE(mesh_sanity_check(*mesh_p)) << "First test mesh with problems!";
}

TEST(lf_edge_create, TestMesh1) {
  using coord_t = Eigen::Vector2d;
  using size_type = mesh::Mesh::size_type;

  // Building the mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

  EXPECT_TRUE(mesh_sanity_check(*mesh_p)) << "Second test mesh with problems!";
}

TEST(lf_hybrid2d, EdgeNumbering) {
  // Construct a one element mesh that consists of a quad:
  MeshFactory mf(2);

  // add nodes
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 0)), 0);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(1, 0)), 1);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(1, 1)), 2);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 1)), 3);

  // add an element
  Eigen::MatrixXd node_coord(2, 4);
  node_coord << 0, 1, 1, 0, 0, 0, 1, 1;
  EXPECT_EQ(0, mf.AddEntity(
                   base::RefEl::kQuad(), {0, 1, 2, 3},
                   std::make_unique<geometry::QuadO1>(std::move(node_coord))));

  // explicitly add the right edge:
  node_coord = Eigen::MatrixXd(2, 2);
  node_coord << 1, 1, 0, 1;
  EXPECT_EQ(0, mf.AddEntity(base::RefEl::kSegment(), {1, 2},
                            std::make_unique<geometry::SegmentO1>(node_coord)));

  // build the mesh
  auto mesh = mf.Build();

  EXPECT_EQ(mesh->NumEntities(0), 1);
  EXPECT_EQ(mesh->NumEntities(1), 4);
  EXPECT_EQ(mesh->NumEntities(2), 4);

  // check indices of the nodes:
  Eigen::VectorXd zero = Eigen::VectorXd::Zero(0);
  auto entities2 = mesh->Entities(2);
  auto node0 = std::find_if(entities2.begin(), entities2.end(), [&](auto& e) {
    return e.Geometry()->Global(zero).norm() < 1e-6;
  });
  EXPECT_NE(node0, entities2.end());
  EXPECT_EQ(mesh->Index(*node0), 0);

  auto node1 = std::find_if(entities2.begin(), entities2.end(), [&](auto& e) {
    return e.Geometry()->Global(zero).isApprox(Eigen::Vector2d(1, 0));
  });
  EXPECT_NE(node1, entities2.end());
  EXPECT_EQ(mesh->Index(*node1), 1);

  auto node2 = std::find_if(entities2.begin(), entities2.end(), [&](auto& e) {
    return e.Geometry()->Global(zero).isApprox(Eigen::Vector2d(1, 1));
  });
  EXPECT_NE(node2, entities2.end());
  EXPECT_EQ(mesh->Index(*node2), 2);

  auto node3 = std::find_if(entities2.begin(), entities2.end(), [&](auto& e) {
    return e.Geometry()->Global(zero).isApprox(Eigen::Vector2d(0, 1));
  });
  EXPECT_NE(node3, entities2.end());
  EXPECT_EQ(mesh->Index(*node3), 3);

  // check index of right edge:
  auto entities1 = mesh->Entities(1);
  auto right_edge =
      std::find_if(entities1.begin(), entities1.end(), [](auto& e) {
        return e.Geometry()
            ->Global(Eigen::VectorXd::Constant(1, 0.5))
            .isApprox(Eigen::Vector2d(1, 0.5));
      });
  EXPECT_NE(right_edge, entities1.end());
  EXPECT_EQ(mesh->Index(*right_edge), 0);

  // check index of element:
  EXPECT_EQ(mesh->Index(*mesh->Entities(0).begin()), 0);
}

TEST(lf_hybrid2d, IncompleteMeshes) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  // construct a mesh with a triangle + an additional point and make sure it
  // fails:
  auto factory = std::make_unique<MeshFactory>(2, true);
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(0, 0)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(1, 0)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(0, 1)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(1, 1)));
  factory->AddEntity(
      base::RefEl::kTria(), {0, 1, 2},
      std::make_unique<geometry::TriaO1>(
          (Eigen::MatrixXd(2, 3) << 0, 1, 0, 0, 0, 1).finished()));

  EXPECT_DEATH(factory->Build(), "Mesh is incomplete");

  // construct the same mesh but with check_completeness = false
  factory = std::make_unique<MeshFactory>(2, false);
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(0, 0)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(1, 0)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(0, 1)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(1, 1)));
  factory->AddEntity(
      base::RefEl::kTria(), {0, 1, 2},
      std::make_unique<geometry::TriaO1>(
          (Eigen::MatrixXd(2, 3) << 0, 1, 0, 0, 0, 1).finished()));
  EXPECT_NO_FATAL_FAILURE(factory->Build());

  // construct a mesh with a triangle + a detached edge:
  factory = std::make_unique<MeshFactory>(2, true);
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(0, 0)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(1, 0)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(0, 1)));
  factory->AddPoint(std::make_unique<geometry::Point>(Eigen::Vector2d(1, 1)));
  factory->AddEntity(
      base::RefEl::kTria(), {0, 1, 2},
      std::make_unique<geometry::TriaO1>(
          (Eigen::MatrixXd(2, 3) << 0, 1, 0, 0, 0, 1).finished()));
  factory->AddEntity(base::RefEl::kSegment(), {2, 3},
                     std::make_unique<geometry::SegmentO1>(
                         (Eigen::Matrix2d() << 0, 1, 1, 1).finished()));
  EXPECT_DEATH(factory->Build(), "Mesh is incomplete");
#pragma GCC diagnostic pop
}

}  // namespace lf::mesh::hybrid2d::test

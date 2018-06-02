#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/geometry/geometry.h>

using namespace lf::mesh;
using namespace lf::geometry;

using size_type = unsigned int;

TEST(hybrid2d, directMeshConstruction) {

  // construct a very simple mesh with two elements (tria, quad)
  // manually.

  std::vector<Eigen::VectorXd> nodes{
    Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 0),
    Eigen::Vector2d(2, 0), Eigen::Vector2d(1, 1), Eigen::Vector2d(0, 1)
  };

  std::vector<std::tuple<std::vector<size_type>, std::unique_ptr<lf::geometry::
                           Geometry>>> elements;
  Eigen::MatrixXd nodesQuad(2,4), nodesTria(2,3);
  nodesQuad << 0, 1, 1, 0,
    0, 0, 1, 1;
  elements.emplace_back(std::vector<size_type>{0, 1, 3, 4},
                        std::make_unique<QuadO1>(nodesQuad));
  nodesTria << 0, 1, 0,
    0, 0, 1;
  elements.emplace_back(std::vector<size_type>{1, 2, 3},
                        std::make_unique<TriaO1>(nodesTria));

  hybrid2d::Mesh mesh(2, std::move(nodes), std::move(elements));

  EXPECT_EQ(mesh.DimMesh(), 2);
  EXPECT_EQ(mesh.DimWorld(), 2);

  EXPECT_EQ(mesh.Size(0), 2);
  EXPECT_EQ(mesh.Size(1), 6);
  EXPECT_EQ(mesh.Size(2), 5);

  auto element0 = mesh.Entities(0).begin();
  auto element1 = element0;
  ++element1;

  auto node0 = mesh.Entities(2).begin();
  auto node1 = node0;
  ++node1;
  auto node2 = node1;
  ++node2;
  auto node3 = node2;
  ++node3;
  auto node4 = node3;
  ++node4;

  EXPECT_EQ(element0->Codim(), 0);
  EXPECT_EQ(element1->Codim(), 0);

  EXPECT_EQ(element0->SubEntities(2)[0], *node0);
  EXPECT_EQ(element0->SubEntities(2)[1], *node1);
  EXPECT_EQ(element0->SubEntities(2)[2], *node3);
  EXPECT_EQ(element0->SubEntities(2)[3], *node4);

  EXPECT_EQ(element0->SubEntities(1)[1], element1->SubEntities(1)[2]);
  EXPECT_EQ(element0->SubEntities(1)[1].SubEntities(1)[0].Codim(), 2);
  EXPECT_EQ(element0->SubEntities(1)[1].SubEntities(1)[0], *node1);
  EXPECT_EQ(element0->SubEntities(1)[1].SubEntities(1)[1], *node3);
}

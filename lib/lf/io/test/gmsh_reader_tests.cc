/**
 * @file
 * @brief Test the GmshReader
 * @author Raffael Casagrande
 * @date   2018-07-01 08:11:24
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/check_entity_indexing.h>
#include <lf/mesh/test_utils/check_geometry_orientation.h>
#include <lf/mesh/test_utils/check_local_topology.h>
#include <lf/mesh/test_utils/check_mesh_completeness.h>

namespace lf::io::test {

using size_type = mesh::Mesh::size_type;

void checkTwoElementMesh(const GmshReader& reader) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  auto mesh = reader.mesh();
  EXPECT_EQ(mesh->NumEntities(0), 2);
  EXPECT_EQ(mesh->NumEntities(1), 6);
  EXPECT_EQ(mesh->NumEntities(2), 5);

  // codim=2 checks:
  auto entities2 = mesh->Entities(2);
  auto origin = std::find_if(entities2.begin(), entities2.end(), [](auto e) {
    return e->Geometry()->Global(Eigen::MatrixXd(0, 1)).squaredNorm() < 1e-10;
  });
  EXPECT_NE(origin, entities2.end());

  EXPECT_EQ(reader.PhysicalEntityNr(**origin).size(), 2);
  EXPECT_EQ(reader.PhysicalEntityNr(**origin)[0], 1);
  EXPECT_EQ(reader.PhysicalEntityNr(**origin)[1], 2);

  EXPECT_EQ(reader.PhysicalEntityNr2Name(1, 2), "physicalEntity1");
  EXPECT_THROW(reader.PhysicalEntityNr2Name(1), base::LfException);
  EXPECT_EQ(reader.PhysicalEntityNr2Name(2, 2), "physicalEntity2");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(2), "physicalEntity2");

  EXPECT_EQ(reader.PhysicalEntityName2Nr("physicalEntity1", 2), 1);
  EXPECT_EQ(reader.PhysicalEntityName2Nr("physicalEntity2", 2), 2);
  EXPECT_THROW(reader.PhysicalEntityName2Nr("physicalEntity1"),
               base::LfException);
  EXPECT_EQ(reader.PhysicalEntityName2Nr("physicalEntity2"), 2);

  EXPECT_THROW(reader.PhysicalEntityNr2Name(100), base::LfException);
  EXPECT_THROW(reader.PhysicalEntityName2Nr("gugus"), base::LfException);

  auto pe2 = reader.PhysicalEntities(2);
  using pe_t = std::pair<size_type, std::string>;
  EXPECT_EQ(pe2.size(), 2);
  EXPECT_NE(std::find(pe2.begin(), pe2.end(), pe_t{1, "physicalEntity1"}),
            pe2.end());
  EXPECT_NE(std::find(pe2.begin(), pe2.end(), pe_t{2, "physicalEntity2"}),
            pe2.end());

  for (auto e : entities2) {
    if (e == *origin) {
      continue;
    }
    EXPECT_EQ(reader.PhysicalEntityNr(*e).size(), 0);
  }

  // codim = 1 checks
  auto entities1 = mesh->Entities(1);
  auto diagonal_edge =
      std::find_if(entities1.begin(), entities1.end(), [](auto e) {
        return e->Geometry()
                   ->Jacobian((Eigen::MatrixXd(1, 1) << 0).finished())
                   .norm() > 1.1;
      });
  EXPECT_NE(diagonal_edge, entities1.end());
  auto diagonal_nr = reader.PhysicalEntityName2Nr("diagonal");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(diagonal_nr), "diagonal");
  EXPECT_EQ(reader.PhysicalEntityNr(**diagonal_edge).size(), 1);
  EXPECT_EQ(reader.PhysicalEntityNr(**diagonal_edge)[0], diagonal_nr);

  auto pe1 = reader.PhysicalEntities(1);
  EXPECT_EQ(pe1.size(), 1);
  EXPECT_EQ(pe1[0].first, 4);
  EXPECT_EQ(pe1[0].second, "diagonal");

  for (auto e : entities1) {
    if (e == *diagonal_edge) {
      continue;
    }
    EXPECT_EQ(reader.PhysicalEntityNr(*e).size(), 0);
  }

  // codim = 0 checks
  auto entities0 = mesh->Entities(0);
  auto square = std::find_if(entities0.begin(), entities0.end(), [](auto e) {
    return e->RefEl() == base::RefEl::kQuad();
  });
  EXPECT_NE(square, entities0.end());
  auto triangle = std::find_if(entities0.begin(), entities0.end(), [](auto& e) {
    return e->RefEl() == base::RefEl::kTria();
  });
  EXPECT_NE(triangle, entities0.end());

  auto square_nr = reader.PhysicalEntityName2Nr("square");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(square_nr), "square");
  EXPECT_EQ(reader.PhysicalEntityNr(**square).size(), 1);
  EXPECT_EQ(reader.PhysicalEntityNr(**square)[0], square_nr);

  EXPECT_EQ(reader.PhysicalEntityName2Nr("physicalEntity1", 0), 1);
  EXPECT_EQ(reader.PhysicalEntityName2Nr("physicalEntity3"), 3);
  EXPECT_EQ(reader.PhysicalEntityNr2Name(1, 0), "physicalEntity1");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(3), "physicalEntity3");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(3, 0), "physicalEntity3");
  EXPECT_THROW(reader.PhysicalEntityNr2Name(3, 1), base::LfException);

  EXPECT_EQ(reader.PhysicalEntityNr(**triangle).size(), 2);
  EXPECT_EQ(reader.PhysicalEntityNr(**triangle)[0], 1);
  EXPECT_EQ(reader.PhysicalEntityNr(**triangle)[1], 3);

  auto pe0 = reader.PhysicalEntities(0);
  EXPECT_EQ(pe0.size(), 3);
  EXPECT_NE(std::find(pe0.begin(), pe0.end(), pe_t{1, "physicalEntity1"}),
            pe0.end());
  EXPECT_NE(std::find(pe0.begin(), pe0.end(), pe_t{3, "physicalEntity3"}),
            pe0.end());
  EXPECT_NE(std::find(pe0.begin(), pe0.end(), pe_t{5, "square"}), pe0.end());

  for (auto e : entities0) {
    mesh::test_utils::checkGeometryOrientation(*e);
    mesh::test_utils::checkLocalTopology(*e);
    mesh::test_utils::checkRelCodim(*e);
  }
  mesh::test_utils::checkEntityIndexing(*reader.mesh());
  mesh::test_utils::checkMeshCompleteness(*reader.mesh());
#pragma GCC diagnostic pop
}

TEST(lf_io, readTwoElementMesh) {
  checkTwoElementMesh(
      GmshReader(std::make_unique<mesh::hybrid2d::MeshFactory>(2),
                 test_utils::getMeshPath("two_element_hybrid_2d.msh")));
  checkTwoElementMesh(
      GmshReader(std::make_unique<mesh::hybrid2d::MeshFactory>(2),
                 test_utils::getMeshPath("two_element_hybrid_2d_binary.msh")));

  checkTwoElementMesh(
      GmshReader(std::make_unique<mesh::hybrid2d::MeshFactory>(2),
                 test_utils::getMeshPath("two_element_hybrid_2d.msh")));

  checkTwoElementMesh(
      GmshReader(std::make_unique<mesh::hybrid2d::MeshFactory>(2),
                 test_utils::getMeshPath("two_element_hybrid_2d_v4.msh")));

  checkTwoElementMesh(GmshReader(
      std::make_unique<mesh::hybrid2d::MeshFactory>(2),
      test_utils::getMeshPath("two_element_hybrid_2d_v4_binary.msh")));

  // Make sure, that we can read a second order mesh:
  auto reader =
      test_utils::getGmshReader("two_element_hybrid_2d_second_order.msh", 2);
  checkTwoElementMesh(reader);
  for (auto e : reader.mesh()->Entities(0)) {
    if (e->RefEl() == base::RefEl::kTria()) {
      EXPECT_TRUE(dynamic_cast<const geometry::TriaO2*>(e->Geometry()));
    } else {
      EXPECT_TRUE(dynamic_cast<const geometry::QuadO2*>(e->Geometry()));
    }
  }

  reader =
      test_utils::getGmshReader("two_element_hybrid_2d_second_order_v4.msh", 2);
  checkTwoElementMesh(reader);
  for (auto e : reader.mesh()->Entities(0)) {
    if (e->RefEl() == base::RefEl::kTria()) {
      EXPECT_TRUE(dynamic_cast<const geometry::TriaO2*>(e->Geometry()));
    } else {
      EXPECT_TRUE(dynamic_cast<const geometry::QuadO2*>(e->Geometry()));
    }
  }
}

TEST(lf_io, readLectureDemoMesh) {
  // the following file contains an extra whitespace at the end of a line...
  // check that we can still read it.
  auto reader = test_utils::getGmshReader("lecturedemomesh.msh", 2);
}

// loads a coarse mesh of a circle parametrized with first order/second order
// elements and computes the volume of the circle and compares it with the exact
// volume.
TEST(lf_io, secondOrderMesh) {
  std::vector<quad::QuadRule> quad_rules(5);
  quad_rules[base::RefEl::kTria().Id()] =
      quad::make_QuadRule(base::RefEl::kTria(), 10);
  quad_rules[base::RefEl::kQuad().Id()] =
      quad::make_QuadRule(base::RefEl::kQuad(), 10);
  auto computeVolume = [&](const mesh::Mesh& m) {
    double result = 0;
    for (auto e : m.Entities(0)) {
      auto& qr = quad_rules[e->RefEl().Id()];
      auto ie = e->Geometry()->IntegrationElement(qr.Points());
      result += ie.cwiseProduct(qr.Weights()).sum();
    }
    return result;
  };

  // with first order triangles:
  auto reader = test_utils::getGmshReader("circle_first_order.msh", 2);
  EXPECT_GT(std::abs(computeVolume(*reader.mesh()) - base::kPi), 0.3);
  reader = test_utils::getGmshReader("circle_first_order_v4.msh", 2);
  EXPECT_GT(std::abs(computeVolume(*reader.mesh()) - base::kPi), 0.3);

  // with triangles
  reader = test_utils::getGmshReader("circle_second_order.msh", 2);
  EXPECT_LT(std::abs(computeVolume(*reader.mesh()) - base::kPi), 0.0032);
  reader = test_utils::getGmshReader("circle_second_order_v4.msh", 2);
  EXPECT_LT(std::abs(computeVolume(*reader.mesh()) - base::kPi), 0.0032);

  // with quadrilaterals:
  reader = test_utils::getGmshReader("circle_second_order_quad.msh", 2);
  EXPECT_LT(std::abs(computeVolume(*reader.mesh()) - base::kPi), 0.003);
  reader = test_utils::getGmshReader("circle_second_order_quad_v4.msh", 2);
  EXPECT_LT(std::abs(computeVolume(*reader.mesh()) - base::kPi), 0.003);
}

void checkPieceOfCake(const GmshReader& reader) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  auto& mesh = *reader.mesh();

  // Total number of entities:
  EXPECT_EQ(mesh.DimMesh(), 2);
  EXPECT_EQ(mesh.DimWorld(), 2);
  EXPECT_EQ(mesh.NumEntities(0), 2);
  EXPECT_EQ(mesh.NumEntities(1), 5);
  EXPECT_EQ(mesh.NumEntities(2), 4);

  // check physical entites:
  auto entities2 = mesh.Entities(2);
  auto origin = std::find_if(entities2.begin(), entities2.end(), [](auto e) {
    return e->Geometry()->Global(Eigen::MatrixXd(0, 1)).squaredNorm() < 1e-10;
  });
  EXPECT_NE(origin, entities2.end());
  EXPECT_EQ(reader.PhysicalEntityNr(**origin), std::vector<unsigned>{1});
  EXPECT_TRUE(reader.IsPhysicalEntity(**origin, 1));

  // get the two curves:
  for (auto e : mesh.Entities(1)) {
    auto nodes = e->Geometry()->Global(base::RefEl::kSegment().NodeCoords());
    if (std::abs(nodes.col(0).norm() - 1) < 1e-6 &&
        std::abs(nodes.col(1).norm() - 1) < 1e-6) {
      EXPECT_EQ(reader.PhysicalEntityNr(*e), std::vector<unsigned>{2});
      EXPECT_TRUE(reader.IsPhysicalEntity(*e, 2));
    } else {
      EXPECT_FALSE(reader.IsPhysicalEntity(*e, 2));
    }
  }

  // make sure the two triangles belong to the physical entity 3:
  for (auto e : mesh.Entities(0)) {
    EXPECT_EQ(reader.PhysicalEntityNr(*e), std::vector<unsigned>{3});
    EXPECT_TRUE(reader.IsPhysicalEntity(*e, 3));
  }

  // mapping to names:
  EXPECT_EQ(reader.PhysicalEntityName2Nr("origin"), 1);
  EXPECT_EQ(reader.PhysicalEntityName2Nr("arc"), 2);
  EXPECT_EQ(reader.PhysicalEntityNr2Name(1), "origin");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(1, 2), "origin");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(2), "arc");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(2, 1), "arc");
  EXPECT_THROW(reader.PhysicalEntityNr2Name(3), base::LfException);
  EXPECT_THROW(reader.PhysicalEntityNr2Name(3, 1), base::LfException);

  // physical entities themselves:
  EXPECT_TRUE(reader.PhysicalEntities(0).empty());
  EXPECT_EQ(reader.PhysicalEntities(1).size(), 1);
  EXPECT_EQ(reader.PhysicalEntities(1)[0].first, 2);
  EXPECT_EQ(reader.PhysicalEntities(1)[0].second, "arc");
  EXPECT_EQ(reader.PhysicalEntities(2).size(), 1);
  EXPECT_EQ(reader.PhysicalEntities(2)[0].first, 1);
  EXPECT_EQ(reader.PhysicalEntities(2)[0].second, "origin");

  for (auto e : mesh.Entities(0)) {
    mesh::test_utils::checkGeometryOrientation(*e);
    mesh::test_utils::checkLocalTopology(*e);
    mesh::test_utils::checkRelCodim(*e);
  }
  mesh::test_utils::checkEntityIndexing(*reader.mesh());
  mesh::test_utils::checkMeshCompleteness(*reader.mesh());

#pragma GCC diagnostic pop
}

TEST(lf_io, pieceOfCake) {
  // make sure when can load a partitioned mesh (with periodic links)
  auto reader = test_utils::getGmshReader("piece_of_cake.msh", 2);
  checkPieceOfCake(reader);
}
}  // namespace lf::io::test

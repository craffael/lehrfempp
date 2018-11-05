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
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/check_entity_indexing.h>
#include <lf/mesh/test_utils/check_geometry_orientation.h>
#include <lf/mesh/test_utils/check_local_topology.h>
#include <lf/mesh/test_utils/check_mesh_completeness.h>

namespace lf::io::test {

using size_type = mesh::Mesh::size_type;

void checkTwoElementMesh(const GmshReader& reader) {
  auto mesh = reader.mesh();
  EXPECT_EQ(mesh->Size(0), 2);
  EXPECT_EQ(mesh->Size(1), 6);
  EXPECT_EQ(mesh->Size(2), 5);

  // codim=2 checks:
  auto entities2 = mesh->Entities(2);
  auto origin = std::find_if(entities2.begin(), entities2.end(), [](auto& e) {
    return e.Geometry()->Global(Eigen::MatrixXd(0, 1)).squaredNorm() < 1e-10;
  });
  EXPECT_NE(origin, entities2.end());

  EXPECT_EQ(reader.PhysicalEntityNr(*origin).size(), 2);
  EXPECT_EQ(reader.PhysicalEntityNr(*origin)[0], 1);
  EXPECT_EQ(reader.PhysicalEntityNr(*origin)[1], 2);

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

  for (auto& e : entities2) {
    if (e == *origin) {
      continue;
    }
    EXPECT_EQ(reader.PhysicalEntityNr(e).size(), 0);
  }

  // codim = 1 checks
  auto entities1 = mesh->Entities(1);
  auto diagonal_edge =
      std::find_if(entities1.begin(), entities1.end(), [](auto& e) {
        return e.Geometry()
                   ->Jacobian((Eigen::MatrixXd(1, 1) << 0).finished())
                   .norm() > 1.1;
      });
  EXPECT_NE(diagonal_edge, entities1.end());
  auto diagonal_nr = reader.PhysicalEntityName2Nr("diagonal");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(diagonal_nr), "diagonal");
  EXPECT_EQ(reader.PhysicalEntityNr(*diagonal_edge).size(), 1);
  EXPECT_EQ(reader.PhysicalEntityNr(*diagonal_edge)[0], diagonal_nr);

  auto pe1 = reader.PhysicalEntities(1);
  EXPECT_EQ(pe1.size(), 1);
  EXPECT_EQ(pe1[0].first, 4);
  EXPECT_EQ(pe1[0].second, "diagonal");

  for (auto& e : entities1) {
    if (e == *diagonal_edge) {
      continue;
    }
    EXPECT_EQ(reader.PhysicalEntityNr(e).size(), 0);
  }

  // codim = 0 checks
  auto entities0 = mesh->Entities(0);
  auto square = std::find_if(entities0.begin(), entities0.end(), [](auto& e) {
    return e.RefEl() == base::RefEl::kQuad();
  });
  EXPECT_NE(square, entities0.end());
  auto triangle = std::find_if(entities0.begin(), entities0.end(), [](auto& e) {
    return e.RefEl() == base::RefEl::kTria();
  });
  EXPECT_NE(triangle, entities0.end());

  auto square_nr = reader.PhysicalEntityName2Nr("square");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(square_nr), "square");
  EXPECT_EQ(reader.PhysicalEntityNr(*square).size(), 1);
  EXPECT_EQ(reader.PhysicalEntityNr(*square)[0], square_nr);

  EXPECT_EQ(reader.PhysicalEntityName2Nr("physicalEntity1", 0), 1);
  EXPECT_EQ(reader.PhysicalEntityName2Nr("physicalEntity3"), 3);
  EXPECT_EQ(reader.PhysicalEntityNr2Name(1, 0), "physicalEntity1");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(3), "physicalEntity3");
  EXPECT_EQ(reader.PhysicalEntityNr2Name(3, 0), "physicalEntity3");
  EXPECT_THROW(reader.PhysicalEntityNr2Name(3, 1), base::LfException);

  EXPECT_EQ(reader.PhysicalEntityNr(*triangle).size(), 2);
  EXPECT_EQ(reader.PhysicalEntityNr(*triangle)[0], 1);
  EXPECT_EQ(reader.PhysicalEntityNr(*triangle)[1], 3);

  auto pe0 = reader.PhysicalEntities(0);
  EXPECT_EQ(pe0.size(), 3);
  EXPECT_NE(std::find(pe0.begin(), pe0.end(), pe_t{1, "physicalEntity1"}),
            pe0.end());
  EXPECT_NE(std::find(pe0.begin(), pe0.end(), pe_t{3, "physicalEntity3"}),
            pe0.end());
  EXPECT_NE(std::find(pe0.begin(), pe0.end(), pe_t{5, "square"}), pe0.end());

  for (auto& e : entities0) {
    mesh::test_utils::checkGeometryOrientation(e);
    mesh::test_utils::checkLocalTopology(e);
  }
  mesh::test_utils::checkEntityIndexing(*reader.mesh());
  mesh::test_utils::checkMeshCompleteness(*reader.mesh());
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
}
}  // namespace lf::io::test

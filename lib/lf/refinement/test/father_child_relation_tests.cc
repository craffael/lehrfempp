/**
 * @file
 * @brief Check that the father child relations that are returned from the
 *        refinement module are correct.
 * @author Raffael Casagrande
 * @date   2018-08-12 01:00:35
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

namespace lf::refinement::test {

void checkFatherChildRelations(const MeshHierarchy& mh,
                               base::size_type father_level) {
  auto father_mesh = mh.getMesh(father_level);
  auto child_mesh = mh.getMesh(father_level + 1);

  // check relations between points:
  auto point_child_infos = mh.PointChildInfos(father_level);
  auto point_father_infos = mh.ParentInfos(father_level + 1, 2);
  auto origin = Eigen::VectorXd::Zero(0);
  auto child_found =
      mesh::utils::make_CodimMeshDataSet<bool>(father_mesh, 2, false);
  for (auto& cp : child_mesh->Entities(2)) {
    auto& father_info = point_father_infos[child_mesh->Index(cp)];
    auto fp = father_info.parent_ptr;
    EXPECT_EQ(father_mesh->Index(*fp), father_info.parent_index);
    EXPECT_EQ(father_info.child_number, 0);
    if (fp->Codim() == 2) {
      // father entity is a point
      EXPECT_TRUE(cp.Geometry()->Global(origin).isApprox(
          fp->Geometry()->Global(origin)));
    } else if (fp->Codim() == 1) {
      // father entity is a segment ->calculate distance to segment:
      auto a = fp->Geometry()->Global(Eigen::VectorXd::Zero(1));
      auto b = fp->Geometry()->Global(Eigen::VectorXd::Ones(1));
      auto c = cp.Geometry()->Global(origin);
      double dist = std::abs((Eigen::MatrixXd(2, 2) << (b - a), (c - a))
                                 .finished()
                                 .determinant()) /
                    (b - a).norm();
      EXPECT_LT(dist, 1e-10);
    }
  }
}

TEST(lf_refinement, FatherChildRelations) {
  auto gmsh_reader =
      io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto base_mesh = gmsh_reader.mesh();

  MeshHierarchy mh(base_mesh,
                   std::make_shared<mesh::hybrid2dp::MeshFactory>(2));

  // mark all edges of the triangle for refinement:
  auto marks = mesh::utils::make_CodimMeshDataSet(base_mesh, 1, false);
  auto triangle = std::find_if(
      base_mesh->Entities(0).begin(), base_mesh->Entities(0).end(),
      [](const auto& e) { return e.RefEl() == base::RefEl::kTria(); });

  for (auto& edge : triangle->SubEntities(1)) {
    (*marks)(edge) = true;
  }

  mh.MarkEdges([&](auto& mesh, auto& e) { return (*marks)(e); });
  mh.RefineMarked();

  {
    // For debug purposes: print the meshs into a vtk file
    io::VtkWriter writer0(base_mesh, "mesh0.vtk");
    io::VtkWriter writer1(mh.getMesh(1), "mesh1.vtk");
  }

  auto child_mesh1 = mh.getMesh(1);
  checkFatherChildRelations(mh, 0);
}

}  // namespace lf::refinement::test

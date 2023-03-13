/**
 * @file
 * @brief Check that the father child relations that are returned from the
 *        refinement module are correct.
 * @author Raffael Casagrande
 * @date   2018-08-12 01:00:35
 * @copyright MIT License
 */

#include "refinement_test_utils.h"

namespace lf::refinement::test {

TEST(lf_refinement, FatherChildRelations) {
  MeshHierarchy::Logger()->set_level(spdlog::level::trace);
  // Read test mesh
  auto gmsh_reader =
      io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto base_mesh = gmsh_reader.mesh();
  // Initialize mesh hierarchy with coarsest mesh
  MeshHierarchy mh(base_mesh, std::make_unique<mesh::hybrid2d::MeshFactory>(2));

  // mark all edges of the triangle for refinement:
  auto marks = mesh::utils::make_CodimMeshDataSet(base_mesh, 1, false);
  const auto triangle = std::find_if(
      base_mesh->Entities(0).begin(), base_mesh->Entities(0).end(),
      [](const auto e) { return e->RefEl() == base::RefEl::kTria(); });
  // Visit all edges of the triangle
  for (const auto *edge : (*triangle)->SubEntities(1)) {
    (*marks)(*edge) = true;
  }
  // Mark all edges of the triangle
  mh.MarkEdges([&](auto & /*mesh*/, auto &e) { return (*marks)(e); });
  // Invoke local refinement
  mh.RefineMarked();
  {
    // For debugging purposes: print the meshs into a vtk file
    io::VtkWriter writer0(base_mesh, "mesh0.vtk");
    io::VtkWriter writer1(mh.getMesh(1), "mesh1.vtk");
  }
  // Obtain pointer to mesh created by local refinement
  auto child_mesh1 = mh.getMesh(1);
  // Checks
  checkFatherChildRelations(mh, 0);
  checkGeometryInParent(mh, 0);
}

}  // namespace lf::refinement::test

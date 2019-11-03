/**
 * @file
 * @brief Doxygen snippets to show lf::assemble::COOMatrix usage
 * @author Ralf Hiptmair
 * @date Tue 29 Oct 2019 03:09:10 PM CET
 * @copyright MIT License
 */

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/refinement/refinement.h>

namespace lf::refinement {
void snippetmh(void) {
  //! [usage]
  // Generate test mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // Construction of a mesh hierarchy requires a factory object
  std::unique_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  // Initialize still flat MESH HIERARCHY containing a single mesh
  lf::refinement::MeshHierarchy multi_mesh(mesh_p, std::move(mesh_factory_ptr));
  // Perform a first step of regular refinement: adds a mesh
  multi_mesh.RefineRegular();
  // For demonstration purposes: Mark edges whose center lies inside a square
  std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &edge)>
      marker =
          [](const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge) -> bool {
    Eigen::MatrixXd ref_c = (Eigen::MatrixXd(1, 1) << 0.5).finished();
    Eigen::VectorXd c(edge.Geometry()->Global(ref_c));
    return ((c[0] > 1.0) && (c[0] < 2.0) && (c[1] > 1.0) && (c[1] < 2.0));
  };
  // mark edges for which predicate returns true
  multi_mesh.MarkEdges(marker);
  // Refine current finest mesh locally: adds another mesh to the hierarchy
  multi_mesh.RefineMarked();
  // Finally obtain pointer to finest mesh in the hierarchy with three meshes
  std::shared_ptr<const mesh::Mesh> fine_mesh = multi_mesh.getMesh(2);
  //! [usage]
}

}  // namespace lf::refinement

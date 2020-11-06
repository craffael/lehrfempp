/**
 * @file
 * @brief Simple demonstrations of mesh refinement
 * @author Ralf Hiptmair
 * @date   March 2019
 * @copyright MIT License
 */

#include "lecturedemorefine.h"

namespace lecturedemo {
using size_type = lf::base::size_type;

// Creation of mesh hierarchy by regular refinement
/* SAM_LISTING_BEGIN_1 */
void regrefMeshSequence(const std::shared_ptr<lf::mesh::Mesh>& mesh_p,
                        int refsteps) {
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, refsteps);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

  // Ouput summary information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes created\n";
  multi_mesh.PrintInfo(std::cout);
  // Number of levels
  size_type L = multi_mesh.NumLevels();

  // Retrieve meshes on all levels
  for (size_type level = 0; level < L; ++level) {
    std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // Reference to current mesh
    const lf::mesh::Mesh& mesh{*lev_mesh_p};
    // Output of mesh information
    std::cout << "==== Mesh on level " << level << ": "
              << mesh.NumEntities(lf::base::RefEl::kPoint()) << " NODEs, "
              << mesh.NumEntities(lf::base::RefEl::kSegment()) << " EDGEs, "
              << mesh.NumEntities(lf::base::RefEl::kTria()) << " TRIAs, "
              << mesh.NumEntities(lf::base::RefEl::kQuad()) << " QUADs."
              << std::endl;
  }
}
/* SAM_LISTING_END_1 */

// Driver function
void lecturedemorefine() {
  std::cout << "LehrFEM++ refinement demos" << std::endl;
  // Obtain hybrid test mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // Request regular refinement, 5 levels
  regrefMeshSequence(mesh_p, 5);
}

}  // namespace lecturedemo

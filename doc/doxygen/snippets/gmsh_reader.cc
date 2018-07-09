/**
 * @file
 * @brief Doxygen snippets to show GmshUsage.
 * @author Raffael Casagrande
 * @date   2018-07-08 07:10:52
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

namespace lf::io {

void foo() {
  //! [usage]
  // Read *.msh_file
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  GmshReader reader(std::move(mesh_factory), "my_mesh.msh");

  // Get Physical Entity Number from it's name (as specified in Gmsh):
  auto air_nr = reader.PhysicalEntityName2Nr("air");

  // Get all codim=0 entities that belong to the "air" physical entity:
  for (auto& e : reader.mesh()->Entities(0)) {
    if (reader.IsPhysicalEntity(e, air_nr)) {
      // This entity belongs to the "air" physical entity.
    }
  }
  //! [usage]
}

}  // namespace lf::io

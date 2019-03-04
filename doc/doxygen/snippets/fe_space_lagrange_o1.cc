/**
 * @file
 * @brief Snippets related to FeSpaceLagrangeO1
 * @author Raffael Casagrande
 * @date   2019-03-04 10:42:54
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::uscalfe {

void foo() {
  //! [usage]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");

  auto fe_space =
      std::make_shared<FeSpaceLagrangeO1<double>>(gmsh_reader.mesh());

  // #dofs = #points in mesh:
  assert(fe_space->LocGlobMap().NoDofs() == gmsh_reader.mesh()->NumEntities(2));
  //! [usage]
}

}  // namespace lf::uscalfe

// This file contains Doxygen snippets for the quick reference dofh document
// It defines a function to hold all code snippets and includes necessary
// imports.

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <iostream>
#include <memory>
#include <span>

void qr_dofh_snippets() {
  {
    //! [dofh DynamicFEDofHandler]
    const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

    auto func = [](const lf::mesh::Entity& e) {
      // return number of dofs associated with e
      return 1;
    };

    lf::assemble::DynamicFEDofHandler dofh(mesh_p, func);
    //! [dofh DynamicFEDofHandler]
  }

  {
    //! [dofh UniformFEDofHandler]
    // Generate a simple test mesh
    const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

    // Initialize a DofHandler for p=2 Lagrange FE
    lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 1},
                                            {lf::base::RefEl::kSegment(), 1},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 1}});
    //! [dofh UniformFEDofHandler]

    //! [dofh Key Methods]
    unsigned num_dofs = dofh.NumDofs();
    //! [dofh Key Methods]

    {
      //! [dofh Key Methods 0]
      const lf::mesh::Entity* e = mesh_p->EntityByIndex(0, 0);

      unsigned num_local_dofs = dofh.NumLocalDofs(*e);
      //! [dofh Key Methods 0]
    }

    {
      //! [dofh Key Methods 1]
      const lf::mesh::Entity* e = mesh_p->EntityByIndex(0, 0);

      std::span<const lf::assemble::gdof_idx_t> idx = dofh.GlobalDofIndices(*e);

      for (auto i : idx) {
        std::cout << i << " ";
      }  // -> prints the global indices of DOFs associated with the entity
      //! [dofh Key Methods 1]
    }

    {
      //! [dofh Key Methods 2]
      const lf::mesh::Entity* e = mesh_p->EntityByIndex(0, 0);

      unsigned num_int_dofs = dofh.NumInteriorDofs(*e);
      //! [dofh Key Methods 2]
    }
  }
}

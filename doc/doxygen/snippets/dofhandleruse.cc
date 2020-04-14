/**
 * @file
 * @brief Doxygen snippets to show usage of LehrFEM++'s geometry functionality
 * @author Ralf Hiptmair
 * @date Jan 28, 2020
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <iostream>

// Code snippet demonstrating the use of DofHandler functionality
// copied from lecturedemodof.cc
//! [pdi]
void printDofInfo(const lf::assemble::DofHandler &dofh) {
  // Obtain pointer to the underlying mesh
  const lf::mesh::Mesh &mesh{*dofh.Mesh()};
  // Number of degrees of freedom managed by the DofHandler object
  const lf::assemble::size_type N_dofs(dofh.NumDofs());
  // Output information about dofs for entities of all co-dimensions
  for (lf::base::dim_t codim = 0; codim <= mesh.DimMesh(); codim++) {
    // Visit all entities of a codimension codim
    for (const lf::mesh::Entity *e : mesh.Entities(codim)) {
      // Fetch unique index of current entity supplied by mesh object
      const lf::base::glb_idx_t e_idx = mesh.Index(*e);
      // Number of shape functions covering current entity
      const lf::assemble::size_type no_dofs(dofh.NumLocalDofs(*e));
      // Obtain global indices of those shape functions ...
      nonstd::span<const lf::assemble::gdof_idx_t> dofarray{
          dofh.GlobalDofIndices(*e)};
      // and print them
      std::cout << *e << ' ' << e_idx << ": " << no_dofs << " dofs = [";
      for (int loc_dof_idx = 0; loc_dof_idx < no_dofs; ++loc_dof_idx) {
        std::cout << dofarray[loc_dof_idx] << ' ';
      }
      std::cout << ']';
      // Also output indices of interior shape functions
      nonstd::span<const lf::assemble::gdof_idx_t> intdofarray{
          dofh.InteriorGlobalDofIndices(*e)};
      std::cout << " int = [";
      for (lf::assemble::gdof_idx_t int_dof : intdofarray) {
        std::cout << int_dof << ' ';
      }
      std::cout << ']' << std::endl;
    }
  }
  // List entities associated with the dofs managed by the current
  // DofHandler object
  for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
    const lf::mesh::Entity &e(dofh.Entity(dof_idx));
    std::cout << "dof " << dof_idx << " -> " << e << " " << mesh.Index(e)
              << std::endl;
  }
}  // end function printDofInfo
//! [pdi]

// Code snippet copied from assembly_tests.cc
//! [ufedofh]
void buildUniformFEDofHandler(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Describe local (uniform!) distribution of degrees of freedom
  // 1 dof per node, 2 dofs per edge, 3 dofs per triangle, 4 dofs per quad
  std::map<lf::base::RefEl, lf::base::size_type> local_dof_distribution{
      {lf::base::RefEl::kPoint(), 1},
      {lf::base::RefEl::kSegment(), 2},
      {lf::base::RefEl::kTria(), 3},
      {lf::base::RefEl::kQuad(), 4}};
  // Construct object encoding local->global dof index mapping
  lf::assemble::UniformFEDofHandler dof_handler(mesh_p, local_dof_distribution);
}
//! [ufedofh]

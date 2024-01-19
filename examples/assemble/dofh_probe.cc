/**
 * @file
 * @brief Probes DofHandler output for one cell of a triangular mesh
 * @author Ralf Hiptmair
 * @date   March 2023
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

#include <iostream>

#include "lf/base/base.h"
#include "lf/base/lf_assert.h"
#include "lf/mesh/entity.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

/**
 * @brief Probes DofHandler output for one cell of a triangular mesh
 */

void dofProbe() {
  std::cout << "Probing DofHandler output for one cell of test mesh #3"
            << '\n';
  // Generate test mesh number 3
  constexpr int selector = 3;
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
  // Create a dof handler object describing a uniform distribution
  // of shape functions
  const lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                         {{lf::base::RefEl::kPoint(), 0},
                                          {lf::base::RefEl::kSegment(), 1},
                                          {lf::base::RefEl::kTria(), 2},
                                          {lf::base::RefEl::kQuad(), 2}});
  // Copious output of information about dof handler
  // PrintInfo(std::cout, dof_handler, 30);
  LF_ASSERT_MSG((mesh_p->NumEntities(0) > 8), "No cell number 8!");
  // Fethc pointer to cell number 8
  const lf::mesh::Entity *cell_p = mesh_p->EntityByIndex(0, 8);
  // Output of DofHandler object
  std::cout << "dofh.NumDofs() = " << dofh.NumDofs() << '\n';
  std::cout << "dofh.NumLocalDofs(*cell_p) = " << dofh.NumLocalDofs(*cell_p)
            << '\n';
  std::cout << "dofh.NumInteriorDofs(*cell_p) = "
            << dofh.NumInteriorDofs(*cell_p) << '\n';
  std::cout << "dofh.GlobalDofIndices(*cell_p) = (";
  const auto gdof_idxs = dofh.GlobalDofIndices(*cell_p);
  for (const lf::base::glb_idx_t idx : gdof_idxs) {
    std::cout << idx << ", ";
  }
  std::cout << ")" << '\n';
  std::cout << "dofh.InteriorGlobalDofIndices(*cell_p) = (";
  const auto idof_idxs = dofh.InteriorGlobalDofIndices(*cell_p);
  for (const lf::base::glb_idx_t idx : idof_idxs) {
    std::cout << idx << ", ";
  }
  std::cout << ")" << '\n';
}

int main(int /*argc*/, char **/*argv*/) {
  dofProbe();
  return 0;
}

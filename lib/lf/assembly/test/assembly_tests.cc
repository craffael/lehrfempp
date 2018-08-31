/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Unit tests for assembly facilities
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <iostream>

#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

#include <lf/assembly/dofhandler.h>

namespace lf::assemble::test {

// Simple test
TEST(lf_assembly, dof_index_test) {
  std::cout << "### TEST: D.o.f. on test mesh" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_EQ(mesh_p->Size(0), 9) << "Test mesh: 9 cells expected!";
  EXPECT_EQ(mesh_p->Size(1), 18) << "Test mesh: 18 edges expected!";
  EXPECT_EQ(mesh_p->Size(2), 10) << "Test mesh: 10 nodes expected!";

  // Describe local (uniform!) distribution of degrees of freedom
  // 1 dof per node, 2 dofs per edge, 3 dofs per triangle, 4 dofs per quad
  LocalStaticDOFs2D local_dof_distribution(1, 2, 3, 4);

  EXPECT_EQ(local_dof_distribution.NoLocDofs(lf::base::RefEl::kPoint()), 1)
      << "1 dof per point expected";
  EXPECT_EQ(local_dof_distribution.NoLocDofs(lf::base::RefEl::kSegment()), 2)
      << "2 dof per edge expected";
  EXPECT_EQ(local_dof_distribution.NoLocDofs(lf::base::RefEl::kTria()), 3)
      << "3 dof per triangle expected";
  EXPECT_EQ(local_dof_distribution.NoLocDofs(lf::base::RefEl::kQuad()), 4)
      << "4 dof per quad expected";
  EXPECT_EQ(local_dof_distribution.TotalNoLocDofs(lf::base::RefEl::kPoint()), 1)
      << "1 shape function per point expected";
  EXPECT_EQ(local_dof_distribution.TotalNoLocDofs(lf::base::RefEl::kSegment()),
            4)
      << "4 shape function  per edge expected";
  EXPECT_EQ(local_dof_distribution.TotalNoLocDofs(lf::base::RefEl::kTria()), 12)
      << "13 shape function  per triangle expected";
  EXPECT_EQ(local_dof_distribution.TotalNoLocDofs(lf::base::RefEl::kQuad()), 16)
      << "16 shape function  per quad expected";

  // Construct dofhandler
  UniformFEDofHandler dof_handler(mesh_p, local_dof_distribution);
  const lf::assemble::size_type N_dofs(dof_handler.GetNoDofs());

  for (lf::base::dim_t codim = 0; codim <= 2; codim++) {
    std::cout << "#### DOFs for entities of co-dimension " << (int)codim
              << std::endl;
    for (const lf::mesh::Entity &e : mesh_p->Entities(codim)) {
      const lf::base::glb_idx_t e_idx = mesh_p->Index(e);
      const lf::assemble::size_type no_dofs(dof_handler.GetNoDofs(e));
      lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t> doflist(
          dof_handler.GetGlobalDofs(e));
      std::cout << e << ' ' << e_idx << ": " << no_dofs << " dofs = [";
      for (const lf::assemble::gdof_idx_t &dof : doflist) {
        EXPECT_LT(dof, N_dofs) << "Dof " << dof << " out of range!";
        std::cout << dof << ' ';
      }
      std::cout << ']' << std::endl;
    }
  }

  std::cout << "#### Entities associated with DOFs" << std::endl;
  for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
    const lf::mesh::Entity &e(dof_handler.GetEntity(dof_idx));
    std::cout << "dof " << dof_idx << " -> " << e << " " << mesh_p->Index(e)
              << std::endl;

    lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t> doflist(
        dof_handler.GetGlobalDofs(e));
    bool dof_found = false;
    for (const lf::assemble::gdof_idx_t &dof : doflist) {
      if (dof == dof_idx) {
        dof_found = true;
        break;
      }
    }
    EXPECT_TRUE(dof_found) << "DOF " << dof_idx << " not in list for " << e
                           << " " << mesh_p->Index(e);
  }

}  // end dof_index_test

}  // namespace lf::assemble::test

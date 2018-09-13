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

#include <lf/assembly/assembler.h>

namespace lf::assemble::test {

/** Rudimentary implementation of an assembler for testing
 *
 * It returns an element matrix with the number of the cell
 * on its diagonal and the negative number of the cell on
 * all other positions.
 */
class TestAssembler {
 public:
  using elem_mat_t = Eigen::Matrix<double, 4, 4>;
  using ElemMat = elem_mat_t &;

  TestAssembler(const lf::mesh::Mesh &mesh) : mesh_(mesh) {}
  bool isActive(const lf::mesh::Entity &) { return true; }
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  elem_mat_t mat_;
  const lf::mesh::Mesh &mesh_;
};

TestAssembler::ElemMat TestAssembler::Eval(const lf::mesh::Entity &cell) {
  const lf::base::glb_idx_t cell_idx = mesh_.Index(cell);
  mat_ = ((double)cell_idx * elem_mat_t::Constant(-1.0));
  mat_.diagonal() = (Eigen::Vector4d::Constant(1.0) * (double)cell_idx);
  return mat_;
}

/** Vector assembler for testing */
class TestVectorAssembler {
 public:
  using elem_vec_t = Eigen::Matrix<double, 4, 1>;
  using ElemVec = elem_vec_t &;

  TestVectorAssembler(const lf::mesh::Mesh &mesh) : mesh_(mesh) {}
  bool isActive(const lf::mesh::Entity &) { return true; }
  ElemVec Eval(const lf::mesh::Entity &cell);

 private:
  elem_vec_t vec_;
  const lf::mesh::Mesh &mesh_;
};

TestVectorAssembler::ElemVec TestVectorAssembler::Eval(
    const lf::mesh::Entity &cell) {
  const lf::base::glb_idx_t cell_idx = mesh_.Index(cell);
  vec_ = ((double)cell_idx * elem_vec_t::Constant(1.0));
  return vec_;
}

// Simple output function for degrees of freedom
void output_dofs_test(const lf::mesh::Mesh &mesh,
                      const lf::assemble::DofHandler &dof_handler) {
  EXPECT_EQ(mesh.DimMesh(), 2) << "Only for 2D meshes";

  const lf::assemble::size_type N_dofs(dof_handler.GetNoDofs());
  for (lf::base::dim_t codim = 0; codim <= 2; codim++) {
    std::cout << "#### DOFs for entities of co-dimension " << (int)codim
              << std::endl;
    for (const lf::mesh::Entity &e : mesh.Entities(codim)) {
      const lf::base::glb_idx_t e_idx = mesh.Index(e);
      const lf::assemble::size_type no_dofs(dof_handler.GetNoLocalDofs(e));
      lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t> doflist(
          dof_handler.GlobalDofIndices(e));
      std::cout << e << ' ' << e_idx << ": " << no_dofs << " dofs = [";
      for (const lf::assemble::gdof_idx_t &dof : doflist) {
        EXPECT_LT(dof, N_dofs) << "Dof " << dof << " out of range!";
        std::cout << dof << ' ';
      }
      std::cout << ']' << std::endl;
    }
  }
}

void output_entities_dofs(const lf::mesh::Mesh &mesh,
                          const lf::assemble::DofHandler &dof_handler) {
  EXPECT_EQ(mesh.DimMesh(), 2) << "Only for 2D meshes";

  const lf::assemble::size_type N_dofs(dof_handler.GetNoDofs());
  std::cout << "#### Entities associated with DOFs" << std::endl;
  for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
    const lf::mesh::Entity &e(dof_handler.GetEntity(dof_idx));
    std::cout << "dof " << dof_idx << " -> " << e << " " << mesh.Index(e)
              << std::endl;
    lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t> doflist(
        dof_handler.GlobalDofIndices(e));
    bool dof_found = false;
    for (const lf::assemble::gdof_idx_t &dof : doflist) {
      if (dof == dof_idx) {
        dof_found = true;
        break;
      }
    }
    EXPECT_TRUE(dof_found) << "DOF " << dof_idx << " not in list for " << e
                           << " " << mesh.Index(e);
  }
}

// Simple test for in initialization of a dof handler
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

  output_dofs_test(*mesh_p, dof_handler);
  output_entities_dofs(*mesh_p, dof_handler);
}  // end dof_index_test

// Same test as above, just with a dynamic dof handler
// Should produce the same output
TEST(lf_assembly, dynamic_dof_index_test) {
  std::cout << "### TEST: Dynamic D.o.f. on test mesh" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_EQ(mesh_p->Size(0), 9) << "Test mesh: 9 cells expected!";
  EXPECT_EQ(mesh_p->Size(1), 18) << "Test mesh: 18 edges expected!";
  EXPECT_EQ(mesh_p->Size(2), 10) << "Test mesh: 10 nodes expected!";

  // Construct dofhandler
  std::function<lf::assemble::size_type(const lf::mesh::Entity &)> dof_layout =
      [](const lf::mesh::Entity &e) -> lf::assemble::size_type {
    switch (e.RefEl()) {
      case lf::base::RefEl::kPoint(): {
        return 1;
      }
      case lf::base::RefEl::kSegment(): {
        return 2;
      }
      case lf::base::RefEl::kTria(): {
        return 3;
      }
      case lf::base::RefEl::kQuad(): {
        return 4;
      }
      default: { LF_ASSERT_MSG(false, "Illegal entity type"); }
    }
  };
  lf::assemble::DynamicFEDofHandler dof_handler(mesh_p, dof_layout);

  output_dofs_test(*mesh_p, dof_handler);
  output_entities_dofs(*mesh_p, dof_handler);
}  // end dof_index_test

// Auxiliary function for testing
void linfe_mat_assembly(const lf::mesh::Mesh &mesh,
                        const lf::assemble::DofHandler &dof_handler) {
  const lf::assemble::size_type N_dofs(dof_handler.GetNoDofs());

  std::cout << N_dofs << " degrees of freedom built" << std::endl;
  EXPECT_EQ(N_dofs, 10) << "Dubious numbers of dofs";

  // Output/store index numbers of entities holding global shape functions
  std::vector<lf::base::glb_idx_t> dof_to_entity_index{};
  for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
    const lf::mesh::Entity &e(dof_handler.GetEntity(dof_idx));
    EXPECT_EQ(e.RefEl(), lf::base::RefEl::kPoint()) << "dofs @ " << e.RefEl();
    const lf::base::glb_idx_t e_idx(mesh.Index(e));
    dof_to_entity_index.push_back(e_idx);
    std::cout << "dof " << dof_idx << " @node " << e_idx << std::endl;
  }

  // Dummy assembler
  TestAssembler assembler{mesh};
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);

  mat = lf::assemble::AssembleMatrixLocally<0, lf::assemble::COOMatrix<double>>(
      dof_handler, assembler);

  std::cout << "Assembled " << mat.rows() << "x" << mat.cols() << " matrix"
            << std::endl;

  // Build sparse matrix from COO format
  Eigen::SparseMatrix<double> Galerkin_matrix(mat.rows(), mat.cols());
  Galerkin_matrix.setFromTriplets(mat.triplets.begin(), mat.triplets.end());

  // Output Galerkin matrix
  std::cout << Galerkin_matrix << std::endl;

  // Convert Galerkin matrix to a dense matrix
  Eigen::MatrixXd dense_Gal_mat = Galerkin_matrix;

  // Reference Galerkin matrix
  Eigen::MatrixXd ref_mat(10, 10);
  ref_mat << 26, -13, -14, -5, 0, 0, -6, -13, -7, -5, -13, 16, -10, -6, -3, 0,
      0, 0, 0, -5, -14, -10, 23, 0, -5, -7, -10, -6, 0, 0, -5, -6, 0, 6, -1, 0,
      0, 0, 0, -5, 0, -3, -5, -1, 6, -3, 0, 0, 0, 0, 0, 0, -7, 0, -3, 7, -4, 0,
      0, 0, -6, 0, -10, 0, 0, -4, 10, -6, 0, 0, -13, 0, -6, 0, 0, 0, -6, 13, -7,
      0, -7, 0, 0, 0, 0, 0, 0, -7, 7, -0, -5, -5, 0, -5, 0, 0, 0, 0, -0, 5;

  for (int i = 0; i < N_dofs; ++i) {
    for (int j = 0; j < N_dofs; ++j) {
      EXPECT_DOUBLE_EQ(dense_Gal_mat(i, j),
                       ref_mat(dof_to_entity_index[i], dof_to_entity_index[j]))
          << " mismatch in entry (" << i << ',' << j << ")";
    }
  }

  // test vector assembler
  TestVectorAssembler vec_assembler{mesh};

  // The actual assembly is done here
  auto rhsvec = lf::assemble::AssembleVectorLocally<0, Eigen::VectorXd>(
      dof_handler, vec_assembler);

  // Output assembled vector
  std::cout << "Assembled vector: " << rhsvec.transpose() << std::endl;

  for (int i = 0; i < N_dofs; ++i) {
    EXPECT_DOUBLE_EQ(rhsvec[i],
                     ref_mat(dof_to_entity_index[i], dof_to_entity_index[i]))
        << " mismatch in component " << i;
  }
}

TEST(lf_assembly, mat_assembly_test) {
  std::cout << "### TEST: Assembly on test mesh" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Construct dofhandler
  // Building a dof handler for linear finite elements
  // Variant I: encoding layout of local dofs in an object
  // lf::assemble::LocalLinearLagrangianFE2D lin_fe_loc_dofs {};
  // lf::assemble::UniformFEDofHandler dof_handler(mesh_p, lin_fe_loc_dofs);
  // Variant II using a map for local dof layout
  lf::assemble::UniformFEDofHandler dof_handler(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});

  DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler << std::endl;

  linfe_mat_assembly(*mesh_p, dof_handler);
}

TEST(lf_assembly, dynamic_dof_test) {
  // Same as the previous test, just based on another version of the dof handler

  std::cout << "### TEST: Assembly based on dynamic dofs" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  std::function<lf::assemble::size_type(const lf::mesh::Entity &)>
      lin_fe_dof_layout =
          [](const lf::mesh::Entity &e) -> lf::assemble::size_type {
    return (e.RefEl().Dimension() == 0) ? 1 : 0;
  };

  lf::assemble::DynamicFEDofHandler dof_handler(mesh_p, lin_fe_dof_layout);

  DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler << std::endl;

  linfe_mat_assembly(*mesh_p, dof_handler);
}

// Dummy assembler adapted to two dof per edge
class EdgeDofAssembler {
 public:
  using elem_mat_t = Eigen::Matrix<double, 8, 8>;
  using ElemMat = elem_mat_t &;

  EdgeDofAssembler(const lf::mesh::Mesh &mesh) : mesh_(mesh) {}
  bool isActive(const lf::mesh::Entity &) { return true; }
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  elem_mat_t mat_;
  const lf::mesh::Mesh &mesh_;
};

EdgeDofAssembler::ElemMat EdgeDofAssembler::Eval(const lf::mesh::Entity &cell) {
  const lf::base::glb_idx_t cell_idx = mesh_.Index(cell);
  mat_ = ((double)cell_idx * elem_mat_t::Constant(-1.0));

  for (int i = 0; i < 8; i += 2) {
    for (int j = 0; j < 8; j += 2) {
      mat_(i, j) = (double)cell_idx;
    }
  }
  mat_.diagonal() =
      (Eigen::Matrix<double, 8, 1>::Constant(1.0) * (double)cell_idx);
  return mat_;
}

// Auxiliary function testing assignment of dofs to edges
void edge_dof_assembly_test(const lf::mesh::Mesh &mesh,
                            const lf::assemble::DofHandler &dof_handler) {
  // Output dofs per entity
  output_dofs_test(mesh, dof_handler);
  output_entities_dofs(mesh, dof_handler);

  const lf::assemble::size_type N_dofs(dof_handler.GetNoDofs());

  std::cout << N_dofs << " degrees of freedom built" << std::endl;
  EXPECT_EQ(N_dofs, 36) << "Dubious numbers of dofs";

  // Output/store index numbers of entities holding global shape functions
  std::vector<lf::base::glb_idx_t> dof_to_entity_index{};
  for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
    const lf::mesh::Entity &e(dof_handler.GetEntity(dof_idx));
    EXPECT_EQ(e.RefEl(), lf::base::RefEl::kSegment()) << "dofs @ " << e.RefEl();
    const lf::base::glb_idx_t e_idx(mesh.Index(e));
    dof_to_entity_index.push_back(e_idx);
    std::cout << "dof " << dof_idx << " @edge " << e_idx << std::endl;
  }

  // Dummy assembler
  EdgeDofAssembler assembler{mesh};
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);

  mat = lf::assemble::AssembleMatrixLocally<0, lf::assemble::COOMatrix<double>>(
      dof_handler, assembler);

  std::cout << "Assembled " << mat.rows() << "x" << mat.cols() << " matrix"
            << std::endl;

  // Build sparse matrix from COO format
  Eigen::SparseMatrix<double> Galerkin_matrix(mat.rows(), mat.cols());
  Galerkin_matrix.setFromTriplets(mat.triplets.begin(), mat.triplets.end());

  // Convert Galerkin matrix to a dense matrix
  Eigen::MatrixXd dense_Gal_mat = Galerkin_matrix;

  // Output Galerkin matrix
  std::cout << dense_Gal_mat << std::endl;

  // Reference Galerkin matrix
  Eigen::MatrixXd ref_mat(36, 36);
  ref_mat << 13, -13, -8, -8, 0, 0, 0, 0, -5, 5, -8, -8, -5, 5, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -13, 13, -8, 8, 0, 0,
      0, 0, -5, -5, -8, 8, -5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5, -5, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, -8, -8, 14, -14, 6, -6, 0, 0, 0, 0, -8, -8, 0, 0,
      0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, -8, 8,
      -14, 14, -6, -6, 0, 0, 0, 0, -8, 8, 0, 0, 0, 0, 0, 0, 0, 0, -6, -6, 0, 0,
      0, 0, 0, 0, 0, 0, -6, -6, 0, 0, 0, 0, 0, 0, 6, -6, 13, -13, -7, -7, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, -7,
      -7, 0, 0, 0, 0, -6, -6, -13, 13, -7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, -6, -6, 0, 0, 0, 0, 0, 0, 0, 0, -6, -6, 7, -7, 0, 0, 0, 0, 0, 0, -7,
      -7, 7, -7, 0, -0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, -7, -7, 0, -0, 0, 0, 0, 0, -7, 7, -7, 7, -0, -0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, -7, -0, -0, -5, -5,
      0, 0, 0, 0, 0, -0, 5, -5, 0, 0, -5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5,
      -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0, 5, -5, 0, 0, 0, 0, -0, -0, -5, 5, 0, 0,
      -5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, -5, 0, 0, 0, 0, 0, 0, 0, 0, -0,
      -0, -8, -8, -8, -8, 0, 0, 0, 0, 0, 0, 10, -10, 0, 0, -2, 2, 2, -2, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8, 8, -8, 8, 0, 0, 0, 0,
      0, 0, -10, 10, 0, 0, -2, -2, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, -5, -5, 0, 0, 0, 0, 0, 0, -5, -5, 0, 0, 6, -6, 1, -1, 0,
      0, 0, 0, 0, 0, 1, -1, -5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, -5, 0, 0,
      0, 0, 0, 0, -5, 5, 0, 0, -6, 6, -1, -1, 0, 0, 0, 0, 0, 0, -1, -1, 5, -5,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, 1, -1,
      3, -3, -2, -2, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, -1, -1, -3, 3, 2, -2, 0, 0, 0, 0, -1,
      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
      -2, 0, 0, -2, 2, 5, -5, -3, -3, 0, 0, 0, 0, 0, 0, -3, -3, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, 0, 0, -2, -2, -5, 5, 3, -3,
      0, 0, 0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 7, -7, -4, -4, 0, 0, 0, 0, 3, -3, -4, -4,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, -3,
      -7, 7, 4, -4, 0, 0, 0, 0, -3, -3, 4, -4, 0, 0, 0, 0, 0, 0, 0, 0, -6, -6,
      -6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 10, -10, 0, 0, 0, 0, 0,
      0, 4, -4, -6, -6, 0, 0, 0, 0, 0, 0, 6, -6, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, -4, -4, -10, 10, 0, 0, 0, 0, 0, 0, -4, -4, 6, -6, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0, 1, -1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 5, -5, 0, 0, 0, 0, 0, 0, -5, 5, 0, 0, -5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5, -5, 0, 0, 0, 0, 0, 0, -5, -5,
      0, 0, -5, -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5, 5, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 3, -3, 0, 0,
      0, 0, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, -3, -3, -3, -3, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 4,
      -4, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -4, -4, -4, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0,
      0, 0, 0, 0, 0, 0, 6, -6, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      -6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, -6, -6, -6, -6, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, -6, 0, 0, 0, 0, 0, 0, 0, 0, -6,
      6, 0, 0, 0, 0, 0, 0, 0, 0, -7, 7, -7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, -7, 0, 0, 0, 0, 0, 0, -7, -7,
      -7, -7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, -7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0, 0, -0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0, 0, 0, 0, 0, 0, 0,
      -0, -0, -0, -0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, -0, 0;

  for (int i = 0; i < N_dofs; ++i) {
    for (int j = 0; j < N_dofs; ++j) {
      EXPECT_DOUBLE_EQ(dense_Gal_mat(i, j), ref_mat(i, j))
          << " mismatch in entry (" << i << ',' << j << ")";
    }
  }
}

TEST(lf_assembly, edge_dof_static) {
  std::cout << "Assembly test for edge dof" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Build a dof handler assigning two dofs to each edge
  lf::assemble::UniformFEDofHandler dof_handler(
      mesh_p, {{lf::base::RefEl::kSegment(), 2}});

  DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler << std::endl;

  edge_dof_assembly_test(*mesh_p, dof_handler);
}

TEST(lf_assembly, edge_dof_dynamic) {
  std::cout << "Assembly test for edge dof" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Build a dof handler assigning two dofs to each edge
  lf::assemble::DynamicFEDofHandler dof_handler(
      mesh_p, [](const lf::mesh::Entity &e) -> lf::assemble::size_type {
        return (e.RefEl() == lf::base::RefEl::kSegment()) ? 2 : 0;
      });

  DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler << std::endl;

  edge_dof_assembly_test(*mesh_p, dof_handler);
}

// Auxiliary assembler class for asembly along the boundary
class BoundaryAssembler {
 public:
  using elem_mat_t = Eigen::Matrix<double, 2, 2>;
  using ElemMat = elem_mat_t &;

  BoundaryAssembler(std::shared_ptr<lf::mesh::Mesh> mesh_p);
  bool isActive(const lf::mesh::Entity &);
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  elem_mat_t mat_;
  std::shared_ptr<lf::mesh::Mesh> mesh_p_;
  lf::mesh::utils::CodimMeshDataSet<lf::base::size_type> cells_at_edges_;
};

BoundaryAssembler::BoundaryAssembler(std::shared_ptr<lf::mesh::Mesh> mesh_p)
    : mesh_p_(std::move(mesh_p)),
      cells_at_edges_(lf::mesh::utils::countNoSuperEntities(mesh_p_, 1, 1)) {}

bool BoundaryAssembler::isActive(const lf::mesh::Entity &edge) {
  LF_VERIFY_MSG(edge.Codim() == 1, "Argument must be edge");
  return (cells_at_edges_(edge) == 1);
}

BoundaryAssembler::ElemMat BoundaryAssembler::Eval(
    const lf::mesh::Entity &edge) {
  const lf::base::glb_idx_t edge_idx = mesh_p_->Index(edge);
  mat_(0, 0) = (double)edge_idx;
  mat_(1, 1) = (double)edge_idx;
  mat_(0, 1) = -(double)edge_idx;
  mat_(1, 0) = -(double)edge_idx;
  return mat_;
}

// Auxiliary function for testing
void linfe_boundary_assembly(std::shared_ptr<lf::mesh::Mesh> mesh_p,
                             const lf::assemble::DofHandler &dof_handler) {
  const lf::assemble::size_type N_dofs(dof_handler.GetNoDofs());

  std::cout << N_dofs << " degrees of freedom in linear FE space" << std::endl;
  EXPECT_EQ(N_dofs, 10) << "Dubious numbers of dofs";

  // Output/store index numbers of entities holding global shape functions
  std::vector<lf::base::glb_idx_t> dof_to_entity_index{};
  for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
    const lf::mesh::Entity &e(dof_handler.GetEntity(dof_idx));
    EXPECT_EQ(e.RefEl(), lf::base::RefEl::kPoint()) << "dofs @ " << e.RefEl();
    const lf::base::glb_idx_t e_idx(mesh_p->Index(e));
    dof_to_entity_index.push_back(e_idx);
    std::cout << "dof " << dof_idx << " @node " << e_idx << std::endl;
  }

  // Test assembler
  BoundaryAssembler assembler{mesh_p};
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);

  // Assembly over edges
  mat = lf::assemble::AssembleMatrixLocally<1, lf::assemble::COOMatrix<double>>(
      dof_handler, assembler);

  std::cout << "Assembled " << mat.rows() << "x" << mat.cols() << " matrix"
            << std::endl;

  // Build sparse matrix from COO format
  Eigen::SparseMatrix<double> Galerkin_matrix(mat.rows(), mat.cols());
  Galerkin_matrix.setFromTriplets(mat.triplets.begin(), mat.triplets.end());

  // Convert Galerkin matrix to a dense matrix
  Eigen::MatrixXd dense_Gal_mat = Galerkin_matrix;

  // Output Galerkin matrix
  std::cout << dense_Gal_mat << std::endl;

  // Reference Galerkin matrix
  Eigen::MatrixXd ref_mat(10, 10);
  ref_mat << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 23, -11, 0, 0, 0, 0, -12, 0, 0, 0, -11,
      24, -13, 0, 0, 0, 0, 0, 0, 0, 0, -13, 27, -14, 0, 0, 0, 0, 0, 0, 0, 0,
      -14, 29, -15, 0, 0, 0, 0, 0, 0, 0, 0, -15, 31, -16, 0, 0, 0, 0, 0, 0, 0,
      0, -16, 33, -17, 0, 0, 0, -12, 0, 0, 0, 0, -17, 29;

  for (int i = 0; i < N_dofs; ++i) {
    for (int j = 0; j < N_dofs; ++j) {
      EXPECT_DOUBLE_EQ(dense_Gal_mat(i, j),
                       ref_mat(dof_to_entity_index[i], dof_to_entity_index[j]))
          << " mismatch in entry (" << i << ',' << j << ")";
    }
  }
}

TEST(lf_assembly, boundary_assembly) {
  std::cout << "### TEST: Assembly along boundary of test mesh" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Construct dofhandler for p.w. linear Lagrangian FE space
  lf::assemble::UniformFEDofHandler dof_handler(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});

  linfe_boundary_assembly(mesh_p, dof_handler);
}

}  // namespace lf::assemble::test

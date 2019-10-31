/**
 * @file
 * @brief Tests for the Product Dofhandler class
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <iostream>

#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

#include <lf/assemble/assemble.h>

#include "../discontinuous_fe_constant.h"
#include "../lagr_fe_quadratic.h"
#include "../loc_comp_dpg.h"
#include "../product_dofhandler.h"
#include "../product_element_matrix_provider_factory.h"
#include "../product_fe_space_factory.h"

#include "assembly_test_utils.h"
namespace projects::dpg::test {

// simple test to check consistencies in the dofhandler:
TEST(product_dof_handler, consistency_test) {
  std::cout << "### Test: component->local numbering:" << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_EQ(mesh_p->NumEntities(0), 9) << "Test mesh: 9 cells expected!";
  EXPECT_EQ(mesh_p->NumEntities(1), 18) << "Test mesh: 18 edges expected!";
  EXPECT_EQ(mesh_p->NumEntities(2), 10) << "Test mesh: 10 nodes expected!";

  // extract entities of all reference types:
  const lf::mesh::Entity& point = *mesh_p->EntityByIndex(2, 0);
  const lf::mesh::Entity& edge = *mesh_p->EntityByIndex(1, 0);
  const lf::mesh::Entity& tria = *mesh_p->EntityByIndex(0, 0);
  const lf::mesh::Entity& quad = *mesh_p->EntityByIndex(0, 5);

  std::cout << "generated mesh" << std::endl;
  // check reference elements:
  EXPECT_EQ(point.RefEl(), lf::base::RefEl::kPoint())
      << "Reference element not point";
  EXPECT_EQ(edge.RefEl(), lf::base::RefEl::kSegment())
      << "Reference element not segment";
  EXPECT_EQ(tria.RefEl(), lf::base::RefEl::kTria())
      << "Reference element not tira";
  EXPECT_EQ(quad.RefEl(), lf::base::RefEl::kQuad())
      << "Reference elment not quad";

  // only one component:
  {
    ProductUniformFEDofHandler dof_handler1(mesh_p,
                                            {{{lf::base::RefEl::kPoint(), 1},
                                              {lf::base::RefEl::kSegment(), 2},
                                              {lf::base::RefEl::kTria(), 3},
                                              {lf::base::RefEl::kQuad(), 4}}});

    std::cout << "constructed first dof handler" << std::endl;
    EXPECT_EQ(dof_handler1.NoInteriorDofs(point), 1);
    EXPECT_EQ(dof_handler1.NoInteriorDofs(edge), 2);
    EXPECT_EQ(dof_handler1.NoInteriorDofs(tria), 3);
    EXPECT_EQ(dof_handler1.NoInteriorDofs(quad), 4);

    EXPECT_EQ(dof_handler1.NoLocalDofs(point), 1)
        << "should be 1, is " << dof_handler1.NoLocalDofs(point);
    EXPECT_EQ(dof_handler1.NoLocalDofs(edge), 4)
        << "should be 1, is " << dof_handler1.NoLocalDofs(edge);
    EXPECT_EQ(dof_handler1.NoLocalDofs(tria), 12);
    EXPECT_EQ(dof_handler1.NoLocalDofs(quad), 16);

    EXPECT_EQ(dof_handler1.NoInteriorDofs(point, 0), 1);
    EXPECT_EQ(dof_handler1.NoInteriorDofs(edge, 0), 2);
    EXPECT_EQ(dof_handler1.NoInteriorDofs(tria, 0), 3);
    EXPECT_EQ(dof_handler1.NoInteriorDofs(quad, 0), 4);

    EXPECT_EQ(dof_handler1.NoLocalDofs(point, 0), 1);
    EXPECT_EQ(dof_handler1.NoLocalDofs(edge, 0), 4);
    EXPECT_EQ(dof_handler1.NoLocalDofs(tria, 0), 12);
    EXPECT_EQ(dof_handler1.NoLocalDofs(quad, 0), 16);
  }

  // two components:
  {
    ProductUniformFEDofHandler dof_handler2(mesh_p,
                                            {{{lf::base::RefEl::kPoint(), 1},
                                              {lf::base::RefEl::kSegment(), 1},
                                              {lf::base::RefEl::kTria(), 2},
                                              {lf::base::RefEl::kQuad(), 2}},

                                             {{lf::base::RefEl::kPoint(), 0},
                                              {lf::base::RefEl::kSegment(), 1},
                                              {lf::base::RefEl::kTria(), 1},
                                              {lf::base::RefEl::kQuad(), 2}}});
    std::cout << "constructed second dof handler" << std::endl;

    // check number of dofs for complete fe space:
    EXPECT_EQ(dof_handler2.NoInteriorDofs(point), 1);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(edge), 2);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(tria), 3);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(quad), 4);

    EXPECT_EQ(dof_handler2.NoLocalDofs(point), 1);
    EXPECT_EQ(dof_handler2.NoLocalDofs(edge), 4);
    EXPECT_EQ(dof_handler2.NoLocalDofs(tria), 12);
    EXPECT_EQ(dof_handler2.NoLocalDofs(quad), 16);

    // check number of dofs for first component:
    EXPECT_EQ(dof_handler2.NoInteriorDofs(point, 0), 1);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(edge, 0), 1);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(tria, 0), 2);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(quad, 0), 2);

    EXPECT_EQ(dof_handler2.NoLocalDofs(point, 0), 1);
    EXPECT_EQ(dof_handler2.NoLocalDofs(edge, 0), 3);
    EXPECT_EQ(dof_handler2.NoLocalDofs(tria, 0), 8);
    EXPECT_EQ(dof_handler2.NoLocalDofs(quad, 0), 10);

    // check number of dofs for second component:
    EXPECT_EQ(dof_handler2.NoInteriorDofs(point, 1), 0);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(edge, 1), 1);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(tria, 1), 1);
    EXPECT_EQ(dof_handler2.NoInteriorDofs(quad, 1), 2);

    EXPECT_EQ(dof_handler2.NoLocalDofs(point, 1), 0);
    EXPECT_EQ(dof_handler2.NoLocalDofs(edge, 1), 1);
    EXPECT_EQ(dof_handler2.NoLocalDofs(tria, 1), 4);
    EXPECT_EQ(dof_handler2.NoLocalDofs(quad, 1), 6);
  }
}

// Simple test for in initialization of a dof handler
TEST(product_dof_handler, dof_index_test) {
  std::cout << "### TEST: D.o.f on test mesh" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_EQ(mesh_p->NumEntities(0), 9) << "Test mesh: 9 cells expected!";
  EXPECT_EQ(mesh_p->NumEntities(1), 18) << "Test mesh: 18 edges expected!";
  EXPECT_EQ(mesh_p->NumEntities(2), 10) << "Test mesh: 10 nodes expected!";

  // Describe local (uniform!) distribution of degrees of freedom
  // 1 dof per node, 2 dofs per edge, 3 dofs per triangle, 4 dofs per quad

  // interpretation as one-component space
  std::vector<std::map<lf::base::RefEl, lf::base::size_type>> local_dof_distr1{
      {{lf::base::RefEl::kPoint(), 1},
       {lf::base::RefEl::kSegment(), 2},
       {lf::base::RefEl::kTria(), 3},
       {lf::base::RefEl::kQuad(), 4}}};

  // Construct dofhandler
  std::cout << "Construction as one component dof handler \n";
  ProductUniformFEDofHandler dof_handler1(mesh_p, local_dof_distr1);
  lf::assemble::test::output_dofs_test(*mesh_p, dof_handler1);
  lf::assemble::test::output_entities_dofs(*mesh_p, dof_handler1);

  // interpretation as two-component space:
  std::vector<std::map<lf::base::RefEl, lf::base::size_type>> local_dof_distr2{
      {{lf::base::RefEl::kSegment(), 1},
       {lf::base::RefEl::kTria(), 0},
       {lf::base::RefEl::kQuad(), 3}},

      {{lf::base::RefEl::kPoint(), 1},
       {lf::base::RefEl::kSegment(), 1},
       {lf::base::RefEl::kTria(), 3},
       {lf::base::RefEl::kQuad(), 1}}};

  ProductUniformFEDofHandler dof_handler2(mesh_p, local_dof_distr2);
  lf::assemble::test::output_dofs_test(*mesh_p, dof_handler2);
  lf::assemble::test::output_entities_dofs(*mesh_p, dof_handler2);

}  // end dof_index_test

TEST(product_dof_handler, mat_assembly_test) {
  std::cout << "### TEST: Assembly on test mesh" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Construct dofhandler using a map for local dof layout
  ProductUniformFEDofHandler dof_handler(mesh_p,
                                         {{{lf::base::RefEl::kPoint(), 1}}});

  lf::assemble::DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler << std::endl;
  lf::assemble::test::linfe_mat_assembly(*mesh_p, dof_handler);
}

TEST(product_dof_handler, edge_dof_static_test) {
  std::cout << "Assembly test for edge dof" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Build a dof handler assigning two dofs to each edge
  ProductUniformFEDofHandler dof_handler1(mesh_p,
                                          {{{lf::base::RefEl::kSegment(), 2}}});

  lf::assemble::DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler1 << std::endl;
  lf::assemble::test::edge_dof_assembly_test(*mesh_p, dof_handler1);
}

TEST(product_dof_handler, boundary_assembly) {
  std::cout << "### TEST: Assembly along boundary of test mesh" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Construct dofhandler for p.w. linear Lagrangian FE space
  ProductUniformFEDofHandler dof_handler(mesh_p,
                                         {{{lf::base::RefEl::kPoint(), 1}}});
  lf::assemble::test::linfe_boundary_assembly(mesh_p, dof_handler);
}

TEST(product_dof_handler, mat_vec_mult_test) {
  std::cout << "### TEST: Local and global multiplication with vector"
            << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Construct dofhandler using a map for local dof layout
  ProductUniformFEDofHandler dof_handler(mesh_p,
                                         {{{lf::base::RefEl::kPoint(), 1}}});

  lf::assemble::DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler << std::endl;
  std::cout << " s= "
            << lf::assemble::test::test_vec_lr_mult(*mesh_p, dof_handler)
            << std::endl;
}

TEST(product_dof_handler, two_component_fe_space) {
  std::cout << "product dofhandler with two zero-th order components. \n "
            << "the expected Galerkin should be block-diagonal. \n"
            << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // construct the dofhandler:
  ProductUniformFESpaceFactory<double> space_fac(mesh_p);
  auto v = space_fac.AddH1Component(1);
  auto u = space_fac.AddL2Component(0);
  auto fe_space = space_fac.Build();

  auto gamma = lf::uscalfe::MeshFunctionConstant(1.0);

  ProductElementMatrixProviderFactory stiffness_factory(fe_space, fe_space);
  stiffness_factory.AddReactionElementMatrixProvider(v, v, gamma);
  stiffness_factory.AddReactionElementMatrixProvider(u, u, gamma);
  auto provider = stiffness_factory.Build();

  lf::assemble::COOMatrix<double> mat(fe_space->LocGlobMap().NoDofs(),
                                      fe_space->LocGlobMap().NoDofs());
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, fe_space->LocGlobMap(), *provider);

  // Build sparse matrix from COO format
  Eigen::SparseMatrix<double> Galerkin_matrix = mat.makeSparse();

  // Convert Galerkin matrix to a dense matrix
  Eigen::MatrixXd dense_Gal_mat = Galerkin_matrix;

  // Output Galerkin matrix
  std::cout << dense_Gal_mat << std::endl;
}

}  // namespace projects::dpg::test

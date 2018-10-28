/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Unit tests for parametric finite element facilities
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "lf/fe/lagr_fe.h"
#include <gtest/gtest.h>
#include <iostream>
#include "lf/fe/loc_comp_ellbvp.h"

#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::fe::test {

TEST(lf_fe, lf_fe_linfe) {
  // Three points in the reference element
  Eigen::MatrixXd refcoords{
      (Eigen::MatrixXd(2, 3) << 0.3, 0.1, 0.7, 0.2, 0.5, 0.1).finished()};
  std::cout << "Points in reference cell\n" << refcoords << std::endl;

  // Testing triangular element
  {
    TriaLinearLagrangeFE<double> tlfe{};
    EXPECT_EQ(tlfe.NumRefShapeFunctions(), 3);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(0, 0), 0);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(1, 0), 0);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(2, 0), 1);

    auto rsf_vals = tlfe.EvalReferenceShapeFunctions(refcoords);
    for (const auto &v : rsf_vals) {
      std::cout << "Tria: RSF values: " << v << std::endl;
    }
    for (const auto &v : tlfe.GradientsReferenceShapeFunctions(refcoords)) {
      std::cout << "Tria: RSF gradients:\n " << v << std::endl;
    }
    std::cout << "Tria: Evaluation nodes\n"
              << tlfe.EvaluationNodes() << std::endl;
  }

  // Testing quadrilateral element
  {
    QuadLinearLagrangeFE<double> qlfe{};
    EXPECT_EQ(qlfe.NumRefShapeFunctions(), 4);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(0, 0), 0);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(1, 0), 0);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(2, 0), 1);

    auto rsf_vals = qlfe.EvalReferenceShapeFunctions(refcoords);
    for (const auto &v : rsf_vals) {
      std::cout << "Quad: RSF values: " << v << std::endl;
    }
    for (const auto &v : qlfe.GradientsReferenceShapeFunctions(refcoords)) {
      std::cout << "Quad: RSF gradients:\n " << v << std::endl;
    }
    std::cout << "Quad: Evaluation nodes\n"
              << qlfe.EvaluationNodes() << std::endl;
  }
}

TEST(lf_fe, lf_fe_ellbvp) {
  std::cout << "### TEST: Computation of element matrices" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  TriaLinearLagrangeFE<double> tlfe{};
  QuadLinearLagrangeFE<double> qlfe{};

  // Set up objects taking care of local computations
  auto alpha = [](Eigen::Vector2d) -> double { return 1.0; };
  auto gamma = [](Eigen::Vector2d) -> double { return 0.0; };
  using loc_comp_t =
      LagrangeFEEllBVPElementMatrix<decltype(alpha), decltype(gamma)>;
  loc_comp_t::ctrl_ = 255;
  lf::quad::QuadRule::out_ctrl_ = 1;

  loc_comp_t comp_elem_mat(tlfe, qlfe, alpha, gamma);

  // For comparison
  LinearFELaplaceElementMatrix lfe_elem_mat{};
  
  // Loop over cells and compute element matrices;
  for (const lf::mesh::Entity &cell : mesh_p->Entities(0)) {
    std::cout << "CELL " << cell << ":" << std::endl;
    std::cout << "Element matrix from LinearFELaplaceElementMatrix:" << std::endl;
    std::cout << lfe_elem_mat.Eval(cell) << std::endl;
    std::cout << "Element matrix from LagrangeFEEllBVPElementMatrix:" << std::endl;
    std::cout << comp_elem_mat.Eval(cell) << std::endl;
  }
}

}  // end namespace lf::fe::test

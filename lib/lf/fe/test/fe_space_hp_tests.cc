/**
 * @file
 * @brief Check that the FeSpaceHierarchic class works as expected
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <cmath>
#include <limits>

namespace lf::fe::test {

TEST(fe_space_hierarchic, legendre_integral) {
  const unsigned N = 1000;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (unsigned p = 0; p < 20; ++p) {
    for (unsigned i = 0; i <= N; ++i) {
      const double x = static_cast<double>(i) / N;
      // Compute the legendre polynomial at x
      const double exact = lf::fe::LegendrePoly<double>::eval(p, x);
      // Differentiate the integral using central differences
      const double intp =
          lf::fe::LegendrePoly<double>::integral(p + 1, x + eps / 2);
      const double intm =
          lf::fe::LegendrePoly<double>::integral(p + 1, x - eps / 2);
      const double approx = (intp - intm) / eps;
      // Compare the two values
      if (std::fabs(exact) < 100) {
        ASSERT_NEAR(approx, exact, 1e-5) << "P=" << p << " x=" << x;
      } else {
        ASSERT_NEAR(approx / exact, 1, 1e-5) << "p=" << p << " x=" << x;
      }
    }
  }
}

TEST(fe_space_hierarchic, legendre_derivative) {
  const unsigned N = 1000;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (unsigned p = 0; p < 20; ++p) {
    for (unsigned i = 0; i <= N; ++i) {
      const double x = static_cast<double>(i) / N;
      // Compute the differentiated legendre polynomial at x
      const double exact = lf::fe::LegendrePoly<double>::derivative(p, x);
      // Differentiate the polynomial using central differences
      const double evalp =
          lf::fe::LegendrePoly<double>::eval(p + 1, x + eps / 2);
      const double evalm =
          lf::fe::LegendrePoly<double>::eval(p + 1, x - eps / 2);
      const double approx = (evalp - evalm) / eps;
      // Compare the two values
      if (std::fabs(exact) < 100) {
        ASSERT_NEAR(approx, exact, 1e-5) << "p=" << p << " x=" << x << " exact=" << exact << " approx=" << approx;
      } else {
        ASSERT_NEAR(approx / exact, 1, 1e-5) << "p=" << p << " x=" << x << " exact=" << exact << " approx=" << approx;
      }
    }
  }
}

TEST(fe_space_hierarchic, jacobi_integral) {
  const unsigned N = 1000;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (unsigned p = 0; p < 20; ++p) {
    for (unsigned alpha = 1; alpha < 10; ++alpha) {
      for (unsigned i = 0; i <= N; ++i) {
        const double x = static_cast<double>(i) / N;
        // Compute the jacobi polynomial at x
        const double exact = lf::fe::JacobiPoly<double>::eval(p, alpha, x);
        // Differentiate the integral using central differences
        const double intp =
            lf::fe::JacobiPoly<double>::integral(p + 1, alpha, x + eps / 2);
        const double intm =
            lf::fe::JacobiPoly<double>::integral(p + 1, alpha, x - eps / 2);
        const double approx = (intp - intm) / eps;
        // Compare the two values
        if (std::fabs(exact) < 100) {
          ASSERT_NEAR(approx, exact, 1e-5)
              << "P=" << p << " alpha=" << alpha << " x=" << x;
        } else {
          ASSERT_NEAR(approx / exact, 1, 1e-5)
              << "p=" << p << " alpha=" << alpha << " x=" << x;
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, continuity) {
  for (int selector = 0; selector <= 8; ++selector) {
    // Get a hybrid test mesh on [0, 3]^2
    const auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    for (unsigned p = 1; p <= 20; ++p) {
      // Generate a FeSpaceHierarchic with order p on the mesh
      const auto fe_space =
          std::make_unique<lf::fe::FeSpaceHierarchic<double>>(mesh, p);
      // Test the continuity of each basis function associated with an edge
      for (auto&& cell : mesh->Entities(0)) {
        const auto sfl_cell = fe_space->ShapeFunctionLayout(*cell);
        const auto cell_refel = cell->RefEl();
        const auto cell_num_nodes = cell_refel.NumNodes();
        const auto cell_nodes = cell_refel.NodeCoords();
        const auto edges = cell->SubEntities(1);
        const auto orient = cell->RelativeOrientations();
        for (int i = 0; i < edges.size(); ++i) {
          const Eigen::RowVectorXd edge_eval_coords =
              Eigen::RowVectorXd::LinSpaced(p + 1, 0, 1);
          const Eigen::MatrixXd cell_eval_coords = [&]() {
            Eigen::MatrixXd result(2, p + 1);
            const Eigen::Vector2d node0 =
                cell_nodes.col((i + 0) % cell_num_nodes);
            const Eigen::Vector2d node1 =
                cell_nodes.col((i + 1) % cell_num_nodes);
            result.row(0) = Eigen::RowVectorXd::Constant(p + 1, node0[0]) +
                            (node1[0] - node0[0]) * edge_eval_coords;
            result.row(1) = Eigen::RowVectorXd::Constant(p + 1, node0[1]) +
                            (node1[1] - node0[1]) * edge_eval_coords;
            return result;
          }();
          const auto edge = edges[i];
          const auto sfl_edge = fe_space->ShapeFunctionLayout(*edge);
          const auto rsf_edge =
              sfl_edge->EvalReferenceShapeFunctions(edge_eval_coords);
          const auto rsf_cell =
              sfl_cell->EvalReferenceShapeFunctions(cell_eval_coords);
          // Compare each basis function on the edge
          ASSERT_TRUE(sfl_edge->NumRefShapeFunctions(0) ==
                      sfl_cell->NumRefShapeFunctions(1))
              << "selector=" << selector << " p=" << p
              << " cell=" << mesh->Index(*cell) << " edge=" << i << std::endl;
          for (int rsf_idx = 0; rsf_idx < sfl_edge->NumRefShapeFunctions(0);
               ++rsf_idx) {
            const Eigen::RowVectorXd rsf_edge_eval = rsf_edge.row(2 + rsf_idx);
            Eigen::RowVectorXd rsf_cell_eval;
            if (orient[i] == lf::mesh::Orientation::positive) {
              rsf_cell_eval =
                  rsf_cell.row(cell_num_nodes + (p - 1) * i + rsf_idx);
            } else {
              rsf_cell_eval =
                  rsf_cell.row(cell_num_nodes + (p - 1) * (i + 1) - rsf_idx - 1)
                      .reverse();
            }
            const double max_diff =
                (rsf_edge_eval - rsf_cell_eval).array().abs().maxCoeff();
            ASSERT_TRUE(max_diff < 1e-10)
                << "selector=" << selector << " p=" << p
                << " cell=" << mesh->Index(*cell) << " edge=" << i
                << " rsf_idx=" << rsf_idx << "\nrsf_edge_eval=["
                << rsf_edge_eval << "]\nrsf_cell_eval=[" << rsf_cell_eval << "]"
                << std::endl;
          }
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, grad_segment) {
  // Test the gradients of the shape functions up to degree p=20
  for (unsigned p = 1; p <= 20; ++p) {
    // Test for all possible combinations of orientations
    lf::mesh::Orientation orient[1];
    for (const auto o0 :
         {lf::mesh::Orientation::positive, lf::mesh::Orientation::negative}) {
      orient[0] = o0;
      const lf::fe::FeHierarchicSegment<double> sfl(p, orient);
      const auto eval_nodes = sfl.EvaluationNodes();
      // Compute the exact gradients
      const auto grad_exact = sfl.GradientsReferenceShapeFunctions(eval_nodes);
      // Compute the gradients with finite differences
      const double dx = std::sqrt(std::numeric_limits<double>::epsilon());
      const auto rsfp = sfl.EvalReferenceShapeFunctions(
          eval_nodes + Eigen::RowVectorXd::Constant(eval_nodes.cols(), dx / 2));
      const auto rsfn = sfl.EvalReferenceShapeFunctions(
          eval_nodes - Eigen::RowVectorXd::Constant(eval_nodes.cols(), dx / 2));
      const Eigen::MatrixXd grad_approx = (rsfp - rsfn) / dx;
      // Compare the gradients
      for (int i = 0; i < grad_exact.rows(); ++i) {
        const Eigen::RowVectorXd grad_i_exact = grad_exact.row(i);
        const Eigen::RowVectorXd grad_i_approx = grad_approx.row(i);
        const double max_diff =
            (grad_i_exact - grad_i_approx).array().abs().maxCoeff();
        ASSERT_TRUE(max_diff < 1e-8) << "i=" << i << std::endl;
      }
    }
  }
}

TEST(fe_space_hierarchic, grad_tria) {
  // Test the gradients of the shape functions up to degree p=20
  for (unsigned p = 1; p <= 20; ++p) {
    // Test for all possible combinations of orientations
    lf::mesh::Orientation orient[4];
    for (const auto o0 :
         {lf::mesh::Orientation::positive, lf::mesh::Orientation::negative}) {
      orient[0] = o0;
      for (const auto o1 :
           {lf::mesh::Orientation::positive, lf::mesh::Orientation::negative}) {
        orient[1] = o1;
        for (const auto o2 : {lf::mesh::Orientation::positive,
                              lf::mesh::Orientation::negative}) {
          orient[2] = o2;
          for (const auto o3 : {lf::mesh::Orientation::positive,
                                lf::mesh::Orientation::negative}) {
            orient[3] = o3;
            const lf::fe::FeHierarchicTria<double> sfl(p, orient);
            const auto eval_nodes = sfl.EvaluationNodes();
            // Compute the exact gradients
            const auto grad_exact =
                sfl.GradientsReferenceShapeFunctions(eval_nodes);
            // Compute the gradients with finite differences
            const double d = std::sqrt(std::numeric_limits<double>::epsilon());
            Eigen::MatrixXd eval_nodes_xp = eval_nodes;
            eval_nodes_xp.row(0) +=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            Eigen::MatrixXd eval_nodes_xn = eval_nodes;
            eval_nodes_xn.row(0) -=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            Eigen::MatrixXd eval_nodes_yp = eval_nodes;
            eval_nodes_yp.row(1) +=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            Eigen::MatrixXd eval_nodes_yn = eval_nodes;
            eval_nodes_yn.row(1) -=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            const auto rsfxp = sfl.EvalReferenceShapeFunctions(eval_nodes_xp);
            const auto rsfxn = sfl.EvalReferenceShapeFunctions(eval_nodes_xn);
            const auto rsfyp = sfl.EvalReferenceShapeFunctions(eval_nodes_yp);
            const auto rsfyn = sfl.EvalReferenceShapeFunctions(eval_nodes_yn);
            const Eigen::MatrixXd dx_approx = (rsfxp - rsfxn) / d;
            const Eigen::MatrixXd dy_approx = (rsfyp - rsfyn) / d;
            // Compare the gradients
            for (int i = 0; i < grad_exact.rows(); ++i) {
              for (int j = 0; j < eval_nodes.cols(); ++j) {
                const double dx_ij_exact = grad_exact(i, 2 * j + 0);
                const double dy_ij_exact = grad_exact(i, 2 * j + 1);
                const double dx_ij_approx = dx_approx(i, j);
                const double dy_ij_approx = dy_approx(i, j);
                ASSERT_TRUE(std::fabs(dx_ij_exact - dx_ij_approx) < 1e-5)
                    << "p=" << p << " i=" << i << " eval_node=["
                    << eval_nodes.col(j).transpose() << "] orient=["
                    << static_cast<int>(orient[0]) << ", "
                    << static_cast<int>(orient[1]) << ", "
                    << static_cast<int>(orient[2]) << "]\ngrad_exact=["
                    << dx_ij_exact << ", " << dy_ij_exact << "]\ngrad_approx=["
                    << dx_ij_approx << ", " << dy_ij_approx << "]" << std::endl;
                ASSERT_TRUE(std::fabs(dy_ij_exact - dy_ij_approx) < 1e-5)
                    << "p=" << p << " i=" << i << " eval_node=["
                    << eval_nodes.col(j).transpose() << "] orient=["
                    << static_cast<int>(orient[0]) << ", "
                    << static_cast<int>(orient[1]) << ", "
                    << static_cast<int>(orient[2]) << "]\ngrad_exact=["
                    << dx_ij_exact << ", " << dy_ij_exact << "]\ngrad_approx=["
                    << dx_ij_approx << ", " << dy_ij_approx << "]" << std::endl;
              }
            }
          }
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, grad_quad) {
  // Test the gradients of the shape functions up to degree p=20
  for (unsigned p = 1; p <= 20; ++p) {
    // Test for all possible combinations of orientations
    lf::mesh::Orientation orient[4];
    for (const auto o0 :
         {lf::mesh::Orientation::positive, lf::mesh::Orientation::negative}) {
      orient[0] = o0;
      for (const auto o1 :
           {lf::mesh::Orientation::positive, lf::mesh::Orientation::negative}) {
        orient[1] = o1;
        for (const auto o2 : {lf::mesh::Orientation::positive,
                              lf::mesh::Orientation::negative}) {
          orient[2] = o2;
          for (const auto o3 : {lf::mesh::Orientation::positive,
                                lf::mesh::Orientation::negative}) {
            orient[3] = o3;
            const lf::fe::FeHierarchicQuad<double> sfl(p, orient);
            const auto eval_nodes = sfl.EvaluationNodes();
            // Compute the exact gradients
            const auto grad_exact =
                sfl.GradientsReferenceShapeFunctions(eval_nodes);
            // Compute the gradients with finite differences
            const double d = std::sqrt(std::numeric_limits<double>::epsilon());
            Eigen::MatrixXd eval_nodes_xp = eval_nodes;
            eval_nodes_xp.row(0) +=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            Eigen::MatrixXd eval_nodes_xn = eval_nodes;
            eval_nodes_xn.row(0) -=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            Eigen::MatrixXd eval_nodes_yp = eval_nodes;
            eval_nodes_yp.row(1) +=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            Eigen::MatrixXd eval_nodes_yn = eval_nodes;
            eval_nodes_yn.row(1) -=
                Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
            const auto rsfxp = sfl.EvalReferenceShapeFunctions(eval_nodes_xp);
            const auto rsfxn = sfl.EvalReferenceShapeFunctions(eval_nodes_xn);
            const auto rsfyp = sfl.EvalReferenceShapeFunctions(eval_nodes_yp);
            const auto rsfyn = sfl.EvalReferenceShapeFunctions(eval_nodes_yn);
            const Eigen::MatrixXd dx_approx = (rsfxp - rsfxn) / d;
            const Eigen::MatrixXd dy_approx = (rsfyp - rsfyn) / d;
            // Compare the gradients
            for (int i = 0; i < grad_exact.rows(); ++i) {
              for (int j = 0; j < eval_nodes.cols(); ++j) {
                const double dx_ij_exact = grad_exact(i, 2 * j + 0);
                const double dy_ij_exact = grad_exact(i, 2 * j + 1);
                const double dx_ij_approx = dx_approx(i, j);
                const double dy_ij_approx = dy_approx(i, j);
                ASSERT_TRUE(std::fabs(dx_ij_exact - dx_ij_approx) < 1e-8)
                    << "p=" << p << " i=" << i << " orient=["
                    << static_cast<int>(orient[0]) << ", "
                    << static_cast<int>(orient[1]) << ", "
                    << static_cast<int>(orient[2]) << ", "
                    << static_cast<int>(orient[3]) << "]\ngrad_exact=["
                    << dx_ij_exact << ", " << dy_ij_exact << "]\ngrad_approx=["
                    << dx_ij_approx << ", " << dy_ij_approx << "]" << std::endl;
                ASSERT_TRUE(std::fabs(dy_ij_exact - dy_ij_approx) < 1e-8)
                    << "p=" << p << " i=" << i << " orient=["
                    << static_cast<int>(orient[0]) << ", "
                    << static_cast<int>(orient[1]) << ", "
                    << static_cast<int>(orient[2]) << ", "
                    << static_cast<int>(orient[3]) << "]\ngrad_exact=["
                    << dx_ij_exact << ", " << dy_ij_exact << "]\ngrad_approx=["
                    << dx_ij_approx << ", " << dy_ij_approx << "]" << std::endl;
              }
            }
          }
        }
      }
    }
  }
}

}  // end namespace lf::fe::test

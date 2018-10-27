/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Unit tests for matrix in COO format
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <iostream>

#include <lf/assemble/fix_dof.h>

namespace lf::assemble::test {

/** Seting up a small matrix in COO format and testing matrix x vector product
 */
TEST(lf_assembly, coomatrix_test) {
  std::cout << "Test for matrix in COO format" << std::endl;

  // Create 8x10 matrix in COO format: all zero in the beginning
  COOMatrix<double> M(8, 10);
  // Fill the matrix

  M.AddToEntry(3 - 1, 1 - 1, -8);
  M.AddToEntry(1 - 1, 2 - 1, 1);
  M.AddToEntry(4 - 1, 2 - 1, -7);
  M.AddToEntry(2 - 1, 3 - 1, 2);
  M.AddToEntry(5 - 1, 3 - 1, -6);
  M.AddToEntry(3 - 1, 4 - 1, 3);
  M.AddToEntry(6 - 1, 4 - 1, -5);
  M.AddToEntry(4 - 1, 5 - 1, 4);
  M.AddToEntry(7 - 1, 5 - 1, -4);
  M.AddToEntry(5 - 1, 6 - 1, 5);
  M.AddToEntry(8 - 1, 6 - 1, -3);
  M.AddToEntry(6 - 1, 7 - 1, 6);
  M.AddToEntry(7 - 1, 8 - 1, 7);
  M.AddToEntry(8 - 1, 9 - 1, 8);

  std::cout << "Matrix M = \n" << M.makeDense() << std::endl;

  Eigen::VectorXd x{Eigen::VectorXd::LinSpaced(10, 0.5, 5.0)};
  Eigen::VectorXd y{Eigen::VectorXd::LinSpaced(10, 9.0, 0.0)};
  Eigen::VectorXd result{Eigen::VectorXd::Zero(8)};
  Eigen::VectorXd exact(8);
  exact << 9, 17, -52, -33, -16, -1, 12, 23;
  Eigen::VectorXd alt = M.makeDense() * (x + y);

  std::cout << "x = " << x.transpose() << ",\ny = " << y.transpose()
            << ",\nresult = " << exact.transpose() << std::endl;

  std::cout << "Matrix x Vector product M*(x+y):" << std::endl;
  std::cout << (M.MatVecMult(x + y)).transpose() << std::endl;
  Eigen::VectorXd diff1 = M.MatVecMult(x + y) - exact;
  EXPECT_DOUBLE_EQ(diff1.norm(), 0.0) << "diff1 = " << diff1.transpose();
  EXPECT_DOUBLE_EQ((M.MatVecMult(x + y) - alt).norm(), 0.0)
      << "Mismatch with dense matrix arithmetic";
  M.MatVecMult(x + y, result);
  std::cout << result.transpose() << std::endl;
  Eigen::VectorXd diff2 = result - exact;
  EXPECT_DOUBLE_EQ(diff2.norm(), 0.0) << "diff2 = " << diff2.transpose();
}

/* MATLAB for comparison
   A = gallery('tridiag',10,-1,2,-1);
   b = (0:9)';
   i = [1 2 4 6 7 8 10];
   ic = [3,5,9];
   fixed = [-1,-2,-3];
   x = A(i,i)\(b(i)-A(i,ic)*fixed'); */

TEST(lf_assembly, fix_dof_test) {
  std::cout << "Test for fixing solution components" << std::endl;
  // Size of matrix
  const int N = 10;
  // Create 8x10 matrix in COO format: all zero in the beginning
  COOMatrix<double> A(N, N);
  // Right-hand-side vector
  Eigen::VectorXd b(N);
  for (int k = 0; k < N; k++) {
    if (k > 0) {
      A.AddToEntry(k, k - 1, -1.0);
    }
    if (k < N - 1) {
      A.AddToEntry(k, k + 1, -1.0);
    }
    A.AddToEntry(k, k, 2);
    b[k] = k + 1;
  }
  std::cout << "Matrix A = \n" << A.makeDense() << std::endl;
  std::cout << "rhs = " << b.transpose() << std::endl;

  // Fixed solution components
  fixed_components_t<double> fixed_solution_components{
      {2, -1.0}, {4, -2.0}, {8, -3.0}};

  fix_solution_components_lse(fixed_solution_components, A, b);
  std::cout << "Modified matrix A = \n" << A.makeDense() << std::endl;
  std::cout << "Modified rhs = " << b.transpose() << std::endl;

  // Solve linear system
  // Initialize sparse matrix
  Eigen::SparseMatrix<double> stiffness_matrix(A.makeSparse());
  // Solve linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(stiffness_matrix);
  Eigen::VectorXd x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    std::cout << "solver failed!" << std::endl;
  }
  std::cout << "Solution x = " << x.transpose() << std::endl;
  Eigen::VectorXd exact(N);
  exact << 1, 1, -1, 0.5, -2, 7.75, 11.5, 8.25, -3, 3.5;
  EXPECT_NEAR((x - exact).norm(), 0.0, 1.0E-12) << "Wrong result!";
}

}  // namespace lf::assemble::test

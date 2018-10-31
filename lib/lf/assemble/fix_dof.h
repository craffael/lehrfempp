#ifndef _LF_FIXDOF_H
#define _LF_FIXDOF_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Support function for fixing solutiuon components of a linear system
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <Eigen/Sparse>
#include "coomatrix.h"

namespace lf::assemble {

/**
 * @brief enforce prescribed solution components
 * @sa fix_solution_components_lse()
 *
 * @tparam SCALAR underlying scalar type, e.g. double
 * @tparam RHSVEC generic vector type for right hand side
 * @param fixed_comp_flags boolean vector whose length agrees with the matrix
 * dimension. A value of `true` indicates that the corresponding solution
 * component is prescribed
 * @param fixed_vec prescribed values are stored in the corresponding
 * components of this vector
 * @param mat reference to the _square_ coefficient matrix in COO format
 * @param rhs reference to the right-hand-side vector
 *
 * ### Requirements for type RHSVECTOR and FIXEDVALVEC
 * An object of type VECTOR or RESULTVECTOR must provide a method `Size()`
 * telling the length of the vector and `operator []` for read/write access to
 * vector entries.
 *
 * ### Algorithm
 *
 * For the sake of simplicity let the prescribed solution components
 * be at the bottom of the solution vector and denote them by
 \f$\widehat{\mathbf{x}}\f$
 * Let the coefficient matrix \f$\mathbf{A}\f$ and the right-hand-side vector
 * be split according to
 * \f[
     \mathbf{A} = \left[\begin{array}{cc} \mathbf{A}_{11} & \mathbf{A}_{12} \\
                                    \mathbf{A}_{21} & \mathbf{A}_{22}
                    \end{array}\right]\quad,\quad
    \mathbf{b} = \left[\begin{array}{c} \mathbf{b}_1 \\ \mathbf{b}_2
 \end{array}\right] \f]
 * Then the solution of the system with fixed solution components is
 * given by
 * \f[
      \left[\begin{array}{cc} \mathbf{A}_{11} & 0  \\
        \mathbf{0} & \mathbf{I}
                    \end{array}\right]
      \left[\begin{array}{c} \mathbf{x} \\ \widehat{\mathbf{x}}
 \end{array}\right]
      =
      \left[\begin{array}{c}
      \mathbf{b}_1 - \mathbf{A}_{12}\widehat{\mathbf{x}} \\
       \widehat{\mathbf{x}}
      \end{array}\right]
   \f]
 */
template <typename SCALAR, typename FIXEDVALVEC, typename RHSVECTOR>
void fix_flagged_solution_components(std::vector<bool> &fixed_comp_flags,
                                     const FIXEDVALVEC &fixed_vec,
                                     COOMatrix<SCALAR> &A, RHSVECTOR &b) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  LF_ASSERT_MSG(N == b.size(),
                "Mismatch N = " << N << " <-> b.size() = " << b.size());
  LF_ASSERT_MSG(N == fixed_comp_flags.size(),
                "Mismatch N = " << N << " <-> fixed_comp_flags.size() = "
                                << fixed_comp_flags.size());
  LF_ASSERT_MSG(
      N == fixed_vec.size(),
      "Mismatch N = " << N << " <-> fixed_vec.size() = " << fixed_vec.size());

  // Multiply sparse matrix with the vector of fixed components and subtract
  // result from right hand side
  A.MatVecMult(-1.0, fixed_vec, b);
  // Set vector components of right-hand-side vector for prescribed values
  for (lf::assemble::gdof_idx_t k = 0; k < N; ++k) {
    if (fixed_comp_flags[k]) {
      b[k] = fixed_vec[k];
    }
  }
  // Set rows and columns of the sparse matrix corresponding to the fixed
  // solution components to zero
  A.setZero([&fixed_comp_flags](gdof_idx_t i, gdof_idx_t j) {
    return (fixed_comp_flags[i] || fixed_comp_flags[j]);
  });
  // lf::assemble::COOMatrix<double>::TripletVec::iterator new_last =
  //     std::remove_if(
  //         A.triplets().begin(), A.triplets().end(),
  //         [&fixed_comp_flags](
  //             typename lf::assemble::COOMatrix<double>::Triplet &triplet) {
  //           return (fixed_comp_flags[triplet.row()] ||
  //                   fixed_comp_flags[triplet.col()]);
  //         });
  // // Adjust size of triplet vector
  // A.triplets().erase(new_last, A.triplets().end());
  // Add Unit diagonal entries corrresponding to fixed components
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    if (fixed_comp_flags[dofnum]) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/**
 * @brief Setting unknowns of a sparse linear system of equations
 *        to fixed values
 * @sa fix_flagged_solution_components()
 *
 * ### Algorithm
 *
 * For the sake of simplicity let the prescribed solution components
 * be at the bottom of the solution vector and denote them by
 \f$\widehat{\mathbf{x}}\f$
 * Let the coefficient matrix \f$\mathbf{A}\f$ and the right-hand-side vector
 * be split according to
 * \f[
 \mathbf{A} = \left[\begin{array}{cc} \mathbf{A}_{11} & \mathbf{A}_{12} \\
 \mathbf{A}_{21} & \mathbf{A}_{22}
 \end{array}\right]\quad,\quad
 \mathbf{b} = \left[\begin{array}{c} \mathbf{b}_1 \\ \mathbf{b}_2
 \end{array}\right] \f]
 * Then the solution of the system with fixed solution components is
 * given by
 * \f[
 \left[\begin{array}{cc} \mathbf{A}_{11} & \mathbf{A}_{12}  \\
 \mathbf{0} & \mathbf{I}
 \end{array}\right]
 \left[\begin{array}{c} \mathbf{x} \\ \widehat{\mathbf{x}}
 \end{array}\right]
 =
 \left[\begin{array}{c}
 \mathbf{b}_1 \\
 \widehat{\mathbf{x}}
 \end{array}\right]
 \f]
*/
template <typename SCALAR, typename FIXEDVALVEC, typename RHSVECTOR>
void fix_flagged_solution_comp_alt(std::vector<bool> &fixed_comp_flags,
                                   const FIXEDVALVEC &fixed_vec,
                                   COOMatrix<SCALAR> &A, RHSVECTOR &b) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  LF_ASSERT_MSG(N == b.size(),
                "Mismatch N = " << N << " <-> b.size() = " << b.size());
  LF_ASSERT_MSG(N == fixed_comp_flags.size(),
                "Mismatch N = " << N << " <-> fixed_comp_flags.size() = "
                                << fixed_comp_flags.size());
  LF_ASSERT_MSG(
      N == fixed_vec.size(),
      "Mismatch N = " << N << " <-> fixed_vec.size() = " << fixed_vec.size());

  // Set vector components of right-hand-side vector for prescribed values
  for (lf::assemble::gdof_idx_t k = 0; k < N; ++k) {
    if (fixed_comp_flags[k]) {
      b[k] = fixed_vec[k];
    }
  }
  // Set rows and columns of the sparse matrix corresponding to the fixed
  // solution components to zero
  A.setZero([&fixed_comp_flags](gdof_idx_t i, gdof_idx_t j) {
    return (fixed_comp_flags[i]);
  });
  // lf::assemble::COOMatrix<double>::TripletVec::iterator new_last =
  //     std::remove_if(
  //         A.triplets().begin(), A.triplets().end(),
  //         [&fixed_comp_flags](
  //             typename lf::assemble::COOMatrix<double>::Triplet &triplet) {
  //           return (fixed_comp_flags[triplet.row()]);
  //         });
  // // Adjust size of triplet vector
  // A.triplets().erase(new_last, A.triplets().end());
  // Add Unit diagonal entries corrresponding to fixed components
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    if (fixed_comp_flags[dofnum]) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/**
 * @brief Information about fixed solution components
 */
template <typename SCALAR>
using fixed_components_t =
    std::vector<std::pair<lf::assemble::gdof_idx_t, SCALAR>>;

/**
 * @brief manipulate a square linear system of equations with a coefficient
 * matrix in COO format so that some solution components attain prescribed
 * values.
 *
 * @tparam SCALAR underlying scalar type, e.g. double
 * @tparam RHSVEC generic vector type for right hand side
 * @param mat reference to the _square_ coefficient matrix in COO format
 * @param rhs reference to the right-hand-side vector
 *
 * ### Requirements for type RHSVECTOR
 * An object of type VECTOR or RESULTVECTOR must provide a method `Size()`
 * telling the length of the vector and `operator []` for read/write access to
 * vector entries.
 *
 * @sa fix_flagged_solution_components()
 */
template <typename SCALAR, typename RHSVECTOR>
void fix_solution_components_lse(
    const fixed_components_t<SCALAR> &fixed_components, COOMatrix<SCALAR> &A,
    RHSVECTOR &b) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  LF_ASSERT_MSG(N == b.size(),
                "Mismatch N = " << N << " <-> b.size() = " << b.size());
  // Set up a vector containing the negative fixed solution components and a
  // flag vector, whose entry is `true` if the corresponding
  // vector component is meant to be fixed
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> fixed_vec(b.size());
  fixed_vec.setZero();
  std::vector<bool> fixed_comp_flags(N, false);
  for (typename fixed_components_t<SCALAR>::const_reference idx_val_pair :
       fixed_components) {
    LF_ASSERT_MSG(idx_val_pair.first < N,
                  "Index " << idx_val_pair.first << " >= N = " << N);
    fixed_vec[idx_val_pair.first] += idx_val_pair.second;
    fixed_comp_flags[idx_val_pair.first] = true;
  }
  fix_flagged_solution_components<SCALAR>(fixed_comp_flags, fixed_vec, A, b);
  // fix_flagged_solution_comp_alt<SCALAR>(fixed_comp_flags, fixed_vec, A, b);
}

}  // namespace lf::assemble
#endif

/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Functions for global finite element assembly of matrices and
 *        right hand side vectors
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef _LF_ASSEMBLE_H
#define _LF_ASSEMBLE_H

#include <Eigen/Sparse>
#include "dofhandler.h"

namespace lf::assemble {
/**
 * @brief A temporaery data structure for matrix in COO format.
 *
 * @tparam basic (usually atomic) scalar type for the matrix
 * Default implementation relies on COO/triplet format.
 */
template <typename SCALAR>
struct COOMatrix {
  using Scalar = SCALAR;
  /** Set up zero matrix of a given size */
  COOMatrix(size_type num_rows, size_type num_cols)
      : n(num_rows), m(num_cols) {}

  COOMatrix(const COOMatrix &) = default;
  COOMatrix(COOMatrix &&) = default;
  COOMatrix &operator=(const COOMatrix &) = default;
  COOMatrix &operator=(COOMatrix &&) = default;

  /**
   * @brief Add a value to the specified entry
   * @param i row index
   * @param j column index
   *
   * Adds another index-value triplet.
   */
  void AddToEntry(gdof_idx_t i, gdof_idx_t j, SCALAR increment) {
    triplets.push_back(Eigen::Triplet<SCALAR>(i, j, increment));
  }

  size_type n, m;                               /**< dimensions of matrix */
  std::vector<Eigen::Triplet<SCALAR>> triplets; /**< COO format data */
};

/**
 * @brief Assembly function for standard assembly of finite element matrices
 *
 * @tparam TMPMATRIX a type fitting the concept of COOMatrix
 * @tparam ASSEMBLER a type providing the computation of element matrices
 * @param dof_handler_trial a dof handler object for column space @see
 * DofHandler
 * @param dof_handler_test a dof handler object for row space @see DofHandler
 * @param assembler assembler object for passing all kinds of data
 * @param matrix matrix object to which the assembled matrix will be added
 */
template <typename TMPMATRIX, class ASSEMBLER>
void AssembleMatrixCellwise(const DofHandler &dof_handler_trial,
                            const DofHandler &dof_handler_test,
                            ASSEMBLER &assembler, TMPMATRIX &matrix) {
  // Type of matrix entries, usually either double or complex.
  using scalar_t = typename TMPMATRIX::Scalar;
  // Type for element matrix
  using elem_mat_t = typename ASSEMBLER::ElemMat;
  // Underlying mesh
  const lf::mesh::Mesh &mesh(dof_handler_trial.getMesh());

  LF_ASSERT_MSG(&mesh == &dof_handler_test.getMesh(),
                "Trial and test space must be defined on the same mesh");

  // Central assembly loop over cells
  for (lf::mesh::Entity &cell : mesh.Entities(0)) {
    // Some cells may be skipped
    if (assembler.isActive(cell)) {
      // Size of element matrix
      const size_type nrows_loc = dof_handler_test.GetNoLocalDofs(cell);
      const size_type ncols_loc = dof_handler_trial.GetNoLocalDofs(cell);
      // row indices of for contributions of cells
      lf::base::RandomAccessRange<const gdof_idx_t> row_idx(
          dof_handler_test.GetGlobalDofs(cell));
      // Column indices of for contributions of cells
      lf::base::RandomAccessRange<const gdof_idx_t> col_idx(
          dof_handler_trial.GetGlobalDofs(cell));
      // Request element matrix from assembler object
      const elem_mat_t elem_mat(assembler.eval(cell));
      LF_ASSERT_MSG(elem_mat.rows() >= nrows_loc,
                    "nrows mismatch " << elem_mat.rows() << " <-> " << nrows_loc
                                      << ", cell " << mesh.Index(cell));
      LF_ASSERT_MSG(elem_mat.cols() >= ncols_loc,
                    "ncols mismatch " << elem_mat.cols() << " <-> " << nrows_loc
                                      << ", cell " << mesh.Index(cell));
      // Assembly double loop
      for (int i = 0; i < nrows_loc; i++) {
        for (int j = 0; j < ncols_loc; j++) {
          matrix.AddToEntry(row_idx[i], col_idx[j], elem_mat(i, j));
        }
      }  // end assembly local double loop
    }    // end if(isActive() )
  }      // end main assembly loop
}  // end AssembleMatrixCellwise

template <typename TMPMATRIX, class ASSEMBLER>
TMPMATRIX AssembleMatrixCellwiseY(const DofHandler &dof_handler_trial,
                                  const DofHandler &dof_handler_test,
                                  ASSEMBLER &assembler) {
  TMPMATRIX matrix(dof_handler_test.GetNoDofs(), dof_handler_trial.GetNoDofs());
  AssembleMatrixCellwise(dof_handler_trial, dof_handler_test, assembler,
                         matrix);
  return matrix;
}
template <typename TMPMATRIX, class ASSEMBLER>
TMPMATRIX AssembleMatrixCellwiseX(const DofHandler &dof_handler,
                                  ASSEMBLER &assembler) {
  return AssembleMatrixCellwiseY<TMPMATRIX, ASSEMBLER>(dof_handler, dof_handler,
                                                       assembler);
}

/** Rudimentary implementation */
class LinearFiniteElementAssembler {
 public:
  using elem_mat_t = Eigen::Matrix<double, 4, 4>;
  using ElemMat = elem_mat_t &;
  LinearFiniteElementAssembler(void) : mat_(elem_mat_t::Constant(1.0)) {}
  bool isActive(lf::mesh::Entity &) { return true; }
  ElemMat eval(lf::mesh::Entity &cell) { return mat_; }
  elem_mat_t mat_;
};

}  // namespace lf::assemble

#endif

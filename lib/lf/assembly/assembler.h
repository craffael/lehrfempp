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
 * @tparam SCALAR basic (usually atomic) scalar type for the matrix
 *
 *
 * This class provides a container for a matrix in triplet
 * format, also known as COO format. It essentially manages
 * a vector of Eigen triplets.
 *
 * #### type requirements for template arguments
 *
 * `SCALAR` must be a type that can serve as a scalar type
 * for Eigen::Matrix.
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

  /** @brief return number of rows */
  Eigen::Index rows(void) const { return n; }
  /** @brief return number of column */
  Eigen::Index cols(void) const { return m; }
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
 * @tparam CODIM co-dimension of mesh entities which should be traversed 
 *               in the course of assembly
 * @tparam TMPMATRIX a type fitting the concept of COOMatrix
 * @tparam ASSEMBLER a type providing the computation of element matrices
 * @param dof_handler_trial a dof handler object for column space @see
 * DofHandler
 * @param dof_handler_test a dof handler object for row space @see DofHandler
 * @param assembler assembler object for passing all kinds of data
 * @param matrix matrix object to which the assembled matrix will be added.
 *               The matrix object is not set to zero in the beginning!
 *
 * This method performs cell-oriented assembly controlled by a local-to-global
 * index map ("dof handler").
 *
 * #### type requirements of template arguments
 *
 * - TMPMATRIX is a rudimentary matrix type and must
 * + provide a constructor taking two matrix dimension arguments
 * + have a method `AddtoEntry(i,j,value_to_add)` for adding to a matrix entry
 * + be copyable and assignable
 * A model type is COOMatrix.
 * - ASSEMBLER is a type capable of local assembly of element matrices. It must
 * + have an `Eval()` method returning the element matrix for a cell
 * + supply an `isActive()` method for selecting cells to be taken into account
 * in assembly
 */
template <int CODIM, typename TMPMATRIX, class ASSEMBLER>
void AssembleMatrixLocally(const DofHandler &dof_handler_trial,
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
  for (const lf::mesh::Entity &entity : mesh.Entities(CODIM)) {
    // Some cells may be skipped
    if (assembler.isActive(entity)) {
      // Size of element matrix
      const size_type nrows_loc = dof_handler_test.GetNoLocalDofs(entity);
      const size_type ncols_loc = dof_handler_trial.GetNoLocalDofs(entity);
      // row indices of for contributions of cells
      lf::base::RandomAccessRange<const gdof_idx_t> row_idx(
          dof_handler_test.GlobalDofIndices(entity));
      // Column indices of for contributions of cells
      lf::base::RandomAccessRange<const gdof_idx_t> col_idx(
          dof_handler_trial.GlobalDofIndices(entity));
      // Request local matrix from assembler object. In the case CODIM = 0,
      // when `entity` is a cell, this is the element matrix
      const elem_mat_t elem_mat(assembler.Eval(entity));
      LF_ASSERT_MSG(elem_mat.rows() >= nrows_loc,
                    "nrows mismatch " << elem_mat.rows() << " <-> " << nrows_loc
                                      << ", entity " << mesh.Index(entity));
      LF_ASSERT_MSG(elem_mat.cols() >= ncols_loc,
                    "ncols mismatch " << elem_mat.cols() << " <-> " << nrows_loc
                                      << ", entity " << mesh.Index(entity));
      // Assembly double loop
      for (int i = 0; i < nrows_loc; i++) {
        for (int j = 0; j < ncols_loc; j++) {
          matrix.AddToEntry(row_idx[i], col_idx[j], elem_mat(i, j));
        }
      }  // end assembly local double loop
    }    // end if(isActive() )
  }      // end main assembly loop
}  // end AssembleMatrixLocally

  template <int CODIM,typename TMPMATRIX, class ASSEMBLER>
TMPMATRIX AssembleMatrixLocally(const DofHandler &dof_handler_trial,
                                 const DofHandler &dof_handler_test,
                                 ASSEMBLER &assembler) {
  TMPMATRIX matrix(dof_handler_test.GetNoDofs(), dof_handler_trial.GetNoDofs());
  AssembleMatrixLocally<CODIM,TMPMATRIX, ASSEMBLER>(
      dof_handler_trial, dof_handler_test, assembler, matrix);
  return matrix;
}
  template <int CODIM,typename TMPMATRIX, class ASSEMBLER>
TMPMATRIX AssembleMatrixLocally(const DofHandler &dof_handler,
                                 ASSEMBLER &assembler) {
  return AssembleMatrixLocally<CODIM,TMPMATRIX, ASSEMBLER>(dof_handler, dof_handler,
                                                      assembler);
}

}  // namespace lf::assemble

#endif

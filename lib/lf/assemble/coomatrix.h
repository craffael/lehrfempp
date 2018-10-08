/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Rudimentary triplet based matrix format
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef _LF_COOMATRIX_H
#define _LF_COOMATRIX_H

#include <Eigen/Sparse>
#include "assembly_types.h"

namespace lf::assemble {
/**
 * @brief A temporaery data structure for matrix in COO format.
 *
 * @tparam SCALAR basic scalar type for the matrix
 *
 *
 * This class provides a container for a matrix in triplet
 * format, also known as COO format see
 * [Wikipedia](https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)).
 * It essentially manages a vector of Eigen triplets_.
 *
 * This matrix format allows efficient incremental matrix construction in the
 * context of local assembly.
 *
 * #### type requirements for template arguments
 *
 * `SCALAR` must be a type that can serve as a scalar type
 * for Eigen::Matrix.
 */
template <typename SCALAR>
class COOMatrix {
 public:
  using Scalar = SCALAR;
  /** Set up zero matrix of a given size */
  COOMatrix(size_type num_rows, size_type num_cols)
      : rows_(num_rows), cols_(num_cols) {}

  COOMatrix(const COOMatrix &) = default;
  COOMatrix(COOMatrix &&) noexcept = default;
  COOMatrix &operator=(const COOMatrix &) = default;
  COOMatrix &operator=(COOMatrix &&) noexcept = default;
  ~COOMatrix() = default;

  /** @brief return number of rows */
  Eigen::Index rows() const { return rows_; }
  /** @brief return number of column */
  Eigen::Index cols() const { return cols_; }
  /**
   * @brief Add a value to the specified entry
   * @param i row index
   * @param j column index
   * @param increment
   *
   * Adds another index-value triplet.
   */
  void AddToEntry(gdof_idx_t i, gdof_idx_t j, SCALAR increment) {
    triplets_.push_back(Eigen::Triplet<SCALAR>(i, j, increment));
  }

  /**
   * @brief Create an Eigen::SparseMatrix out of the COO format.
   * @return The created sparse matrix
   *
   * @note This method can be called multiple times and it does not modify the
   * data stored in this COOMatrix.
   */
  Eigen::SparseMatrix<Scalar> makeSparse() const {
    Eigen::SparseMatrix<Scalar> result(rows_, cols_);
    result.setFromTriplets(triplets_.begin(), triplets_.end());
    return result;
  }

 private:
  size_type rows_, cols_;                        /**< dimensions of matrix */
  std::vector<Eigen::Triplet<SCALAR>> triplets_; /**< COO format data */
};

}  // namespace lf::assemble

#endif

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
 * @tparam SCALAR basic (usually atomic) scalar type for the matrix
 *
 *
 * This class provides a container for a matrix in triplet
 * format, also known as COO format see
 * [Wikipedia](https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)).
 * It essentially manages a vector of Eigen triplets.
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

}  // namespace lf::assemble

#endif

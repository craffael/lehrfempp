#ifndef _LF_COOMATRIX_H
#define _LF_COOMATRIX_H
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
  using Triplet = Eigen::Triplet<SCALAR>;
  using TripletVec = std::vector<Triplet>;
  using Index = Eigen::Index;

  /** Set up zero matrix of a given size */
  COOMatrix(size_type num_rows, size_type num_cols)
      : rows_(num_rows), cols_(num_cols) {}

  COOMatrix(const COOMatrix &) = default;
  COOMatrix(COOMatrix &&) noexcept = default;
  COOMatrix &operator=(const COOMatrix &) = default;
  COOMatrix &operator=(COOMatrix &&) noexcept = default;
  ~COOMatrix() = default;

  /** @brief return number of rows */
  Index rows() const { return rows_; }
  /** @brief return number of column */
  Index cols() const { return cols_; }
  /**
   * @brief Add a value to the specified entry
   * @param i row index
   * @param j column index
   * @param increment
   *
   * Adds another index-value triplet.
   *
   * @note the size information for the matrix is adjusted to the passed
   * indices. This ensures that the internally stored numbers of columns and
   * rows are always bigger than the indices of any triplet. When manipulating
   * the triplet information directly, one has to make sure that the size
   * information remains valid.
   */
  void AddToEntry(gdof_idx_t i, gdof_idx_t j, SCALAR increment) {
    rows_ = (i + 1 > rows_) ? i + 1 : rows_;
    cols_ = (j + 1 > cols_) ? j + 1 : cols_;
    triplets_.push_back(Eigen::Triplet<SCALAR>(i, j, increment));
  }
  /**
   * @brief Erase all entries of the matrix
   *
   * This method clears the vector of triplets, effectively setting the
   * matrix to zero. It does not affect the size information about the matrix.
   */
  void setZero() { triplets_.clear(); }
  /**
   * @brief Erase specific entries of the COO matrix, that is, set them to zero
   * @tparam PREDICATE a predicate type compliant with
   *                   std::function<bool(gdof_idx_t,gdof_idx_t)>
   * @param pred selector predicate
   *
   * Removes all triplets `trp` for which `pred(trp.row(),trp.col())` is `true`.
   * This amounts to setting the corresponding matrix entries to zero.
   */
  template <typename PREDICATE>
  void setZero(PREDICATE &&pred) {
    const typename TripletVec::iterator new_last = std::remove_if(
        triplets_.begin(), triplets_.end(),
        [pred](Triplet &trp) { return (pred(trp.row(), trp.col())); });
    // Adjust size of triplet vector
    triplets_.erase(new_last, triplets_.end());
  }

  /**
   * @brief Gives access to the vector of triplets
   * @return reference to internal vector of triplets
   *
   * Use of this method is deprecated. Use setZero(pred) and AddToEntry()
   * instead.
   */
  TripletVec &triplets() { return triplets_; }
  const TripletVec &triplets() const { return triplets_; }

  /**
   * @brief Computes the product of a vector with the matrix in COO format
   * @tparam VECTOR a basic vector type for the argument vector
   * @param alpha scalar with which to multiply the argument vector
   *        before the matrix x vector multiplication.
   * @param vec constant reference to a vector of type VECTOR
   * @return result vector, a dense vector of Eigen
   *
   * ### Requirements for type VECTOR
   * An object of type VECTOR must provide a method `Size()` telling the length
   * of the vector and `operator []` for read access to vector entries.
   */
  template <typename VECTOR>
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> MatVecMult(SCALAR alpha,
                                                      const VECTOR &vec) const;

  /**
   * @brief In-situ computes the product of a vector with the matrix in COO
   * format
   * @tparam VECTOR a basic vector type for the argument vector
   * @tparam RESULTVECTOR another vector type for returning the result
   * @param alpha scalar with which to multiply the argument vector
   *        before the matrix x vector multiplication.
   * @param vec constant reference to a vector of type VECTOR
   * @param result reference to an object of type RESULTVECTOR.
   *
   * @note the result of the matrix-vector product will be added to the entries
   * of `result`!
   *
   * ### Requirements for types VECTOR and RESULTVECTOR
   * An object of type VECTOR or RESULTVECTOR must provide a method `Size()`
   * telling the length of the vector and `operator []` for read/write access to
   * vector entries.
   */
  template <typename VECTOR, typename RESULTVECTOR>
  void MatVecMult(SCALAR alpha, const VECTOR &vec, RESULTVECTOR &resvec) const;

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
  /**
   * @brief Create an Eigen::MatrixX from the COO format
   * @return A dense matrix representing the COO matrix
   *
   * This method is mainly meant for debugging
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> makeDense() const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> mat(rows_, cols_);
    mat.setZero();
    for (const Triplet &trp : triplets_) {
      LF_ASSERT_MSG(
          ((trp.col() < cols_) && (trp.row() < rows_)),
          "Illegal triplet index (" << trp.col() << ',' << trp.row() << ")");
      mat(trp.row(), trp.col()) += trp.value();
    }
    return mat;
  }

  /** @brief Ouput of triplet list describing COO matrix
   *
   * @parm o output stream
   */
  void PrintInfo(std::ostream &o) const {
    o << rows_ << " x " << cols_ << " COO matrix" << std::endl;
    for (const Triplet &trp : triplets_) {
      o << "(" << trp.row() << ',' << trp.col() << ") -> " << trp.value()
        << std::endl;
    }
  }

  /**
   * @brief Output operator for COO matrix
   *
   * This function prints matrix size and the list of triplets
   */
  template <typename SCALARTYPE>
  friend std::ostream &operator<<(std::ostream &o,
                                  const COOMatrix<SCALARTYPE> &mat);

 private:
  size_type rows_, cols_; /**< dimensions of matrix */
  TripletVec triplets_;   /**< COO format data */
};

template <typename SCALARTYPE>
std::ostream &operator<<(std::ostream &o, const COOMatrix<SCALARTYPE> &mat) {
  mat.PrintInfo();
  return o;
}

template <typename SCALAR>
template <typename VECTOR>
Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> COOMatrix<SCALAR>::MatVecMult(
    SCALAR alpha, const VECTOR &vec) const {
  LF_ASSERT_MSG(vec.size() >= cols_,
                "size mismatch: " << cols_ << " <-> " << vec.size());
  Eigen::VectorXd result(rows_);
  result.setZero();
  for (const Triplet &trp : triplets_) {
    result[trp.row()] += trp.value() * (alpha * vec[trp.col()]);
  }
  return result;
}

template <typename SCALAR>
template <typename VECTOR, typename RESULTVECTOR>
void COOMatrix<SCALAR>::MatVecMult(SCALAR alpha, const VECTOR &vec,
                                   RESULTVECTOR &resvec) const {
  LF_ASSERT_MSG(vec.size() >= cols_,
                "Vector vec size mismatch: " << cols_ << " <-> " << vec.size());
  LF_ASSERT_MSG(
      result.size() >= rows_,
      "Vector result size mismatch: " << cols_ << " <-> " << result.size());
  for (const Triplet &trp : triplets_) {
    result[trp.row()] += trp.value() * (alpha * vec[trp.col()]);
  }
}

}  // namespace lf::assemble

#endif

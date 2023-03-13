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

#ifndef INCG_LF_ASSEMBLE_H
#define INCG_LF_ASSEMBLE_H

#include <fmt/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <iostream>

#include "dofhandler.h"

namespace lf::assemble {

/**
 * @defgroup assemble_matrix_locally Cell-Oriented Assembly of Galerkin Matrices
 * @brief Based on helper objects that provide element matrices these functions
 * rely on the distribute assembly policy to set the entries of global Galerkin
 * matrices.
 *
 * All the functions have the name `AssembleMatrixLocally`, but they differ in
 * their arguments. What they have in common is that they are all templated with
 * two types:
 * - ENTITY_MATRIX_PROVIDER is a @ref entity_matrix_provider
 * - TMPMATRIX, a rudimentary matrix type with an `AddToEntry(` mutbale access
 * function
 *
 * The principles of local cell-oriented assembly of Galerkin matrices are
 * explained in [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{sss:assalg}, in particular [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{par:assmat}. An example for the use of local assembly is given in
 * [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{ex:gallse}.
 *
 *  @{
 */

/**
 * @brief The logger that is used by AssembleMatrixLocally() to log additional
 * information. (for logging levels trace + debug)
 */
std::shared_ptr<spdlog::logger> &AssembleMatrixLogger();

/**
 * @brief Assembly function for standard assembly of finite element matrices
 *
 * @tparam TMPMATRIX a type fitting the concept of COOMatrix
 * @tparam ENTITY_MATRIX_PROVIDER a type providing the computation of element
 * matrices, must model the concept \ref entity_matrix_provider
 * @param codim co-dimension of mesh entities which should be traversed
 *              in the course of assembly
 * @param dof_handler_trial a dof handler object for _column space_, see @ref
 * DofHandler
 * @param dof_handler_test a dof handler object for _row space_, see @ref
 * DofHandler
 * @param entity_matrix_provider @ref entity_matrix_provider object for passing
 * all kinds of data
 * @param matrix matrix object to which the assembled matrix will be added.
 *
 * The rationale for passing two different @ref DofHandler objects for trial and
 * test space is given in [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{rem:pg}.
 *
 * @note The matrix object passed in `matrix` is not set to zero in the
 * beginning! This makes is possible to assemble a matrix piecemeal via several
 * successive calls to @ref AssembleMatrixLocally(), see
 * [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{cpp:lfelaplbdmat}.
 *
 * This method performs cell-oriented assembly controlled by a local-to-global
 * index map ("dof handler").
 *
 * #### Logger
 * This function logs additional information to \ref AssembleMatrixLogger(). See
 * \ref loggers for more information.
 *
 * #### type requirements of template arguments
 *
 * - TMPMATRIX is a rudimentary matrix type and must
 * + provide a constructor taking two matrix dimension arguments
 * + have a method `AddtoEntry(i,j,value_to_add)` for adding to a matrix entry
 * A model type is COOMatrix.
 * - ENTITY_MATRIX_PROVIDER is a \ref entity_matrix_provider
 *
 * @note The element matrix returned by the `Eval()` method of
 * `entity_matrix_provider` may have a size larger than that suggested by the
 * number of local shape functions. In this case only its upper left block is
 * accessed.
 *
 * #### Example Usage
 * @snippet assembler.cc matrix_usage
 *
 */
template <typename TMPMATRIX, class ENTITY_MATRIX_PROVIDER>
void AssembleMatrixLocally(dim_t codim, const DofHandler &dof_handler_trial,
                           const DofHandler &dof_handler_test,
                           ENTITY_MATRIX_PROVIDER &entity_matrix_provider,
                           TMPMATRIX &matrix) {
  // Fetch pointer to underlying mesh
  auto mesh = dof_handler_trial.Mesh();
  LF_ASSERT_MSG(mesh == dof_handler_test.Mesh(),
                "Trial and test space must be defined on the same mesh");
  // Central assembly loop over entities of co-dimension specified by
  // the function argument codim
  for (const lf::mesh::Entity *entity : mesh->Entities(codim)) {
    // Some entities may be skipped
    if (entity_matrix_provider.isActive(*entity)) {
      // log the entity reference element + it's global index on level Debug
      SPDLOG_LOGGER_DEBUG(AssembleMatrixLogger(), "Entity {} ({})", *entity,
                          mesh->Index(*entity));
      // Size, aka number of rows and columns, of element matrix
      const size_type nrows_loc = dof_handler_test.NumLocalDofs(*entity);
      const size_type ncols_loc = dof_handler_trial.NumLocalDofs(*entity);
      // row indices of for contributions of cells
      nonstd::span<const gdof_idx_t> row_idx(
          dof_handler_test.GlobalDofIndices(*entity));
      // Column indices of for contributions of cells
      nonstd::span<const gdof_idx_t> col_idx(
          dof_handler_trial.GlobalDofIndices(*entity));
      // Request local matrix from entity_matrix_provider object. In the
      // case codim = 0, when `entity` is a cell, this is the element matrix
      const auto elem_mat{entity_matrix_provider.Eval(*entity)};
      LF_ASSERT_MSG(elem_mat.rows() >= nrows_loc,
                    "nrows mismatch " << elem_mat.rows() << " <-> " << nrows_loc
                                      << ", entity " << mesh->Index(*entity));
      LF_ASSERT_MSG(elem_mat.cols() >= ncols_loc,
                    "ncols mismatch " << elem_mat.cols() << " <-> " << nrows_loc
                                      << ", entity " << mesh->Index(*entity));

      // Log global row and column indices on level debug
      SPDLOG_LOGGER_DEBUG(AssembleMatrixLogger(), "row_idx = {}, col_idx = {}",
                          row_idx, col_idx);

      // Log shape of element matrix on level debug
      SPDLOG_LOGGER_DEBUG(AssembleMatrixLogger(), "{} x {} element matrix",
                          nrows_loc, ncols_loc);

      // Log element matrix itself on level trace
      // TODO(craffael): Doesn't work yet because of fmt issue:
      // https://github.com/fmtlib/fmt/issues/1885
      // SPDLOG_LOGGER_TRACE(AssembleMatrixLogger, elem_mat);

      // Assembly double loop
      std::stringstream ss;  // used to log all triplets on one line.

      for (int i = 0; i < nrows_loc; i++) {
        for (int j = 0; j < ncols_loc; j++) {
          // Add the element at position (i,j) of the local matrix
          // to the entry at (row_idx[i], col_idx[j]) of the global matrix
          matrix.AddToEntry(row_idx[i], col_idx[j], elem_mat(i, j));

          // log the added "triplet" on level trace:
          if (AssembleMatrixLogger()->should_log(spdlog::level::trace)) {
            // if we are on level trace, build string of all triplets:
            ss << "(" << row_idx[i] << ',' << col_idx[j]
               << ")+= " << elem_mat(i, j) << ", ";
          }
        }
      }  // end assembly local double loop

      // log all the triplets on one line on level trace:
      SPDLOG_LOGGER_TRACE(AssembleMatrixLogger(), ss.str());

    }  // end if(isActive() )
  }    // end main assembly loop
}  // end AssembleMatrixLocally

/**
 * @brief Entity-wise local assembly of a matrix from local matrices
 * @tparam TMPMATRIX a type fitting the concept of COOMatrix
 * @tparam ENTITY_MATRIX_PROVIDER a type providing the computation of element
 * matrices, must model the concept \ref entity_matrix_provider
 *
 * @return assembled matrix in a format determined by the template argument
 *         TPMATRIX
 * @sa  AssembleMatrixLocally(const DofHandler &dof_handler_trial,const
 * DofHandler &dof_handler_test,ENTITY_MATRIX_PROVIDER &entity_matrix_provider,
 * TMPMATRIX &matrix)
 *
 * @note Extra requirements for the type TMPMATRIX is imposed; it must
 *       provide the method `setZero()` for setting all entries of the
 *       matrix to zero. It must also possess a constructor that takes
 *       row and column numbers and creates an (empty) matrix of that size.
 *
 * #### Logger
 * This function logs additional information to \ref AssembleMatrixLogger(). See
 * \ref loggers for more information.
 */
template <typename TMPMATRIX, class ENTITY_MATRIX_PROVIDER>
TMPMATRIX AssembleMatrixLocally(
    dim_t codim, const DofHandler &dof_handler_trial,
    const DofHandler &dof_handler_test,
    ENTITY_MATRIX_PROVIDER &entity_matrix_provider) {
  TMPMATRIX matrix{dof_handler_test.NumDofs(), dof_handler_trial.NumDofs()};
  matrix.setZero();
  AssembleMatrixLocally<TMPMATRIX, ENTITY_MATRIX_PROVIDER>(
      codim, dof_handler_trial, dof_handler_test, entity_matrix_provider,
      matrix);
  return matrix;
}

/**
 * @brief Entity-wise local assembly of a matrix from local matrices
 * @tparam TMPMATRIX a type fitting the concept of COOMatrix
 * @tparam ENTITY_MATRIX_PROVIDER a type providing the computation of element
 * matrices, must model the concept \ref entity_matrix_provider
 *
 * @return assembled matrix in a format determined by the template argument
 *         TPMATRIX
 *
 * This special version of the function should be used whenever test and trial
 * space are the same.
 *
 * @sa  AssembleMatrixLocally(const DofHandler &dof_handler_trial,const
 * DofHandler &dof_handler_test,ENTITY_MATRIX_PROVIDER &element_matrix_provider,
 * TMPMATRIX &matrix)
 *
 * #### Logger
 * This function logs additional information to \ref AssembleMatrixLogger(). See
 * \ref loggers for more information.
 */

template <typename TMPMATRIX, class ENTITY_MATRIX_PROVIDER>
TMPMATRIX AssembleMatrixLocally(
    dim_t codim, const DofHandler &dof_handler,
    ENTITY_MATRIX_PROVIDER &entity_matrix_provider) {
  return AssembleMatrixLocally<TMPMATRIX, ENTITY_MATRIX_PROVIDER>(
      codim, dof_handler, dof_handler, entity_matrix_provider);
}
/** @} */  // end of group assemble_matrix_locally

/**
 * @defgroup assemble_vector_locally Cell-Oriented Assembly of R.h.s. Galerkin
 Vectors
 * @brief Functions `assembleVectorLocally()` for setting entries of the vector
 representing the finite element Galerkin discretization of a linear functional.
 *
 * The functions differ in their arguments, but all are templated with the
 following types:
 * - ENTITY_VECTOR_PROVIDER matching the concept @ref entity_vector_provider
 * - VECTOR a generic vector type compliant with Eigen::Vector
 *
 * Also refer to [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{par:vectorassembler} and [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{ex:gallse}.
 *
 * @{
 */
/**
 * @brief entity-local assembly of (right-hand-side) vectors from element
 * vectors
 *
 * @tparam VECTOR a generic vector type with component access through []
 * @tparam ENTITY_VECTOR_PROVIDER type for objects computing entity-local
 * vectors, models concept \ref entity_vector_provider
 * @param codim co-dimension of entities over which assembly should be carried
 * out
 * @param dof_handler object providing local-to-global dof index mapping, see
 * DofHandler
 * @param entity_vector_provider local entity_vector_provider object (passed as
 * non-const!)
 * @param resultvector generic vector for returning the assembled vector
 *
 * ### Type requirements for template arguments
 *
 * - VECTOR must provide a `size()` method telling its length and
 *   read/write access through the `[]` operator.
 *
 * @note Contributions of element vectors are added to the entries of the
 *       `resultvector` argument. This means that `resultvector` has to be
 *       initialized before calling this function!
 *
 * ### Example Usage
 * @snippet assembler.cc vector_usage
 */
template <typename VECTOR, class ENTITY_VECTOR_PROVIDER>
void AssembleVectorLocally(dim_t codim, const DofHandler &dof_handler,
                           ENTITY_VECTOR_PROVIDER &entity_vector_provider,
                           VECTOR &resultvector) {
  // Pointer to underlying mesh
  auto mesh = dof_handler.Mesh();

  // Central assembly loop over entities of the co-dimension specified via
  // the template argument CODIM
  for (const lf::mesh::Entity *entity : mesh->Entities(codim)) {
    // Some cells may be skipped
    if (entity_vector_provider.isActive(*entity)) {
      // Length of element vector
      const size_type veclen = dof_handler.NumLocalDofs(*entity);
      // global dof indices for contribution of the entity
      nonstd::span<const gdof_idx_t> dof_idx(
          dof_handler.GlobalDofIndices(*entity));
      // Request local vector from entity_vector_provider object. In the case
      // CODIM = 0, when `entity` is a cell, this is the element vector
      const auto elem_vec{entity_vector_provider.Eval(*entity)};
      LF_ASSERT_MSG(elem_vec.size() >= veclen,
                    "length mismatch " << elem_vec.size() << " <-> " << veclen
                                       << ", entity " << mesh->Index(*entity));
      // Assembly (single) loop
      for (int i = 0; i < veclen; i++) {
        resultvector[dof_idx[i]] += elem_vec[i];
      }  // end assembly localloop
    }    // end if(isActive() )
  }      // end main assembly loop
}  // end AssembleVectorLocally

/**
 * @brief entity-local assembly of (right-hand-side) vectors from element
 * vectors
 * @tparam VECTOR a generic vector type with component access through []
 * @tparam ENTITY_VECTOR_PROVIDER type for objects computing entity-local
 * vectors, models concept \ref entity_vector_provider
 * @return assembled vector as an object of a type specified by the
 *         VECTOR template argument
 * @sa AssembleVectorLocally(const DofHandler &dof_handler,
 *                           ENTITY_VECTOR_PROVIDER &entity_vector_provider,
 * VECTOR &resultvector)
 *
 * ### Additional type requirements for VECTOR template argument
 *
 * VECTOR must supply a `setZero` method for initialization with zero.
 * All matrix types of Eigen have such a method see
 * [Eigen
 * documentation](https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html)
 */
template <typename VECTOR, class ENTITY_VECTOR_PROVIDER>
VECTOR AssembleVectorLocally(dim_t codim, const DofHandler &dof_handler,
                             ENTITY_VECTOR_PROVIDER &entity_vector_provider) {
  // Allocated vector holding r.h.s. vector to be assembled
  VECTOR resultvector{dof_handler.NumDofs()};
  // Initialize to zero: assembly of new vector
  resultvector.setZero();
  // Perform actual assembly
  AssembleVectorLocally<VECTOR, ENTITY_VECTOR_PROVIDER>(
      codim, dof_handler, entity_vector_provider, resultvector);
  return resultvector;
}  // end AssembleVectorLocally

/** @} */  // end group assemble_vector_locally

}  // namespace lf::assemble

#endif

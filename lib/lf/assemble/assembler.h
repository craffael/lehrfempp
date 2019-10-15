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

#include <iostream>

#include "dofhandler.h"

namespace lf::assemble {

/** @brief debugging output control for matrix assembly function
 *
 * Supported debugging output:
 * - bit 1: prints information about the entity (`amd_entity`)
 * - bit 2: prints information about size of element matrix (`amd_lmdim`)
 * - bit 3: prints element matrix (`md_locmat`)
 * - bit 4: prints global indices of local shape functions (`amd_gdof`)
 */
// extern unsigned int ass_mat_dbg_ctrl;
const unsigned int amd_entity = 1;
const unsigned int amd_lmdim = 2;
const unsigned int amd_locmat = 4;
const unsigned int amd_gdof = 8;
const unsigned int amd_lass = 16;

// The following macro implements:
// extern int ass_mat_dbg_ctrl;
// static lf::base::StaticVar ctrlvarass_mat_dbg_ctrl(
//     "Assembly_ctrl", ass_mat_dbg_ctrl, lf::base::ctrl_root,
//     "Debugging output control for AssembleMatrixLocally()");
// EXTERNDECLAREINFO(ass_mat_dbg_ctrl, "Assembly_ctrl",
//                  "Debugging output control for AssembleMatrixLocally()");

extern unsigned int ass_mat_dbg_ctrl;

/**
 * @brief Assembly function for standard assembly of finite element matrices
 *
 * @tparam TMPMATRIX a type fitting the concept of COOMatrix
 * @tparam ENTITY_MATRIX_PROVIDER a type providing the computation of element
 * matrices
 * @param codim co-dimension of mesh entities which should be traversed
 *              in the course of assembly
 * @param dof_handler_trial a dof handler object for column space @see
 * DofHandler
 * @param dof_handler_test a dof handler object for row space @see DofHandler
 * @param entity_matrix_provider entity_matrix_provider object for passing all
 * kinds of data
 * @param matrix matrix object to which the assembled matrix will be added.
 *
 * @note The matrix object passed in `matrix` is not set to zero in the
 * beginning! This makes is possible to assemble a matrix piecemeal via several
 * successive calls to @ref AssembleMatrixLocally()
 *
 * This method performs cell-oriented assembly controlled by a local-to-global
 * index map ("dof handler").
 *
 * #### type requirements of template arguments
 *
 * - TMPMATRIX is a rudimentary matrix type and must
 * + provide a constructor taking two matrix dimension arguments
 * + have a method `AddtoEntry(i,j,value_to_add)` for adding to a matrix entry
 * A model type is COOMatrix.
 * - ENTITY_MATRIX_PROVIDER is a type capable of local assembly of element
 * matrices. It must
 * + have an `Eval()` method returning the element matrix for a cell
 * + supply an `isActive()` method for selecting cells to be taken into account
 * in assembly
 * - The return type of `Eval()` must be a dense matrix type modelled after
 *   Eigen::Matrix. It must provide:
 * + methods `rows()` and `cols()` telling the number of rows and columns
 * + access to entries via `operator (int,int) const`
 *
 * @note The element matrix returned by the `Eval()` method of
 * `entity_matrix_provider` may have a size larger than that suggested by the
 * number of local shape functions. In this case only its upper left block is
 * accessed.
 */
template <typename TMPMATRIX, class ENTITY_MATRIX_PROVIDER>
void AssembleMatrixLocally(dim_t codim, const DofHandler &dof_handler_trial,
                           const DofHandler &dof_handler_test,
                           ENTITY_MATRIX_PROVIDER &entity_matrix_provider,
                           TMPMATRIX &matrix) {
  // Pointer to underlying mesh
  auto mesh = dof_handler_trial.Mesh();

  LF_ASSERT_MSG(mesh == dof_handler_test.Mesh(),
                "Trial and test space must be defined on the same mesh");

  // Central assembly loop over entities of co-dimension specified by
  // the function argument codim
  for (const lf::mesh::Entity &entity : mesh->Entities(codim)) {
    // Some entities may be skipped
    if (entity_matrix_provider.isActive(entity)) {
      SWITCHEDSTATEMENT(ass_mat_dbg_ctrl, amd_entity,
                        std::cout << "ASM: " << entity << '('
                                  << mesh->Index(entity) << ')' << std::endl);
      // Size, aka number of rows and columns, of element matrix
      const size_type nrows_loc = dof_handler_test.NoLocalDofs(entity);
      const size_type ncols_loc = dof_handler_trial.NoLocalDofs(entity);
      // row indices of for contributions of cells
      nonstd::span<const gdof_idx_t> row_idx(
          dof_handler_test.GlobalDofIndices(entity));
      // Column indices of for contributions of cells
      nonstd::span<const gdof_idx_t> col_idx(
          dof_handler_trial.GlobalDofIndices(entity));
      // Request local matrix from entity_matrix_provider object. In the case
      // codim = 0, when `entity` is a cell, this is the element matrix
      const auto elem_mat{entity_matrix_provider.Eval(entity)};
      LF_ASSERT_MSG(elem_mat.rows() >= nrows_loc,
                    "nrows mismatch " << elem_mat.rows() << " <-> " << nrows_loc
                                      << ", entity " << mesh->Index(entity));
      LF_ASSERT_MSG(elem_mat.cols() >= ncols_loc,
                    "ncols mismatch " << elem_mat.cols() << " <-> " << nrows_loc
                                      << ", entity " << mesh->Index(entity));
      // clang-format off
      SWITCHEDSTATEMENT(
          ass_mat_dbg_ctrl, amd_gdof,
	  std::cout << "ASM: row_idx = ";
          for (auto gdof_idx: row_idx) {
	    std::cout << gdof_idx << ' ';
	  }
	  std::cout << std::endl << "ASM: col_idx = ";
	  for (auto gdof_idx : col_idx) {
            std::cout << gdof_idx << ' ';
          }
	  std::cout << std::endl);
      // clang-format on
      SWITCHEDSTATEMENT(ass_mat_dbg_ctrl, amd_lmdim,
                        std::cout << "ASM: " << nrows_loc << " x " << ncols_loc
                                  << " element matrix" << std::endl);
      SWITCHEDSTATEMENT(
          ass_mat_dbg_ctrl, amd_locmat, for (int i = 0; i < nrows_loc; i++) {
            std::cout << "[ ";
            for (int j = 0; j < ncols_loc; j++) {
              std::cout << elem_mat(i, j) << ' ';
            }
            std::cout << "]" << std::endl;
          });
      // Assembly double loop
      for (int i = 0; i < nrows_loc; i++) {
        for (int j = 0; j < ncols_loc; j++) {
          // Add the element at position (i,j) of the local matrix
          // to the entry at (row_idx[i], col_idx[j]) of the global matrix
          matrix.AddToEntry(row_idx[i], col_idx[j], elem_mat(i, j));
          SWITCHEDSTATEMENT(ass_mat_dbg_ctrl, amd_lass,
                            std::cout << "(" << row_idx[i] << ',' << col_idx[j]
                                      << ")+= " << elem_mat(i, j) << ", ";);
        }
      }  // end assembly local double loop
      SWITCHEDSTATEMENT(ass_mat_dbg_ctrl, amd_lass, std::cout << std::endl;);
    }  // end if(isActive() )
  }    // end main assembly loop
}  // end AssembleMatrixLocally

/**
 * @brief Entity-wise local assembly of a matrix from local matrices
 *
 * @return assembled matrix in a format determined by the template argument
 *         TPMATRIX
 * @sa  AssembleMatrixLocally(const DofHandler &dof_handler_trial,const
 * DofHandler &dof_handler_test,ENTITY_MATRIX_PROVIDER &entity_matrix_provider,
 * TMPMATRIX &matrix)
 *
 * @note An extra requirement for the type TMPMATRIX is imposed; it must
 *       provide the method `setZero()` for setting all entries of the
 *       matrix to zero.
 */
template <typename TMPMATRIX, class ENTITY_MATRIX_PROVIDER>
TMPMATRIX AssembleMatrixLocally(
    dim_t codim, const DofHandler &dof_handler_trial,
    const DofHandler &dof_handler_test,
    ENTITY_MATRIX_PROVIDER &entity_matrix_provider) {
  TMPMATRIX matrix{dof_handler_test.NoDofs(), dof_handler_trial.NoDofs()};
  matrix.setZero();
  AssembleMatrixLocally<TMPMATRIX, ENTITY_MATRIX_PROVIDER>(
      codim, dof_handler_trial, dof_handler_test, entity_matrix_provider,
      matrix);
  return matrix;
}

/**
 * @brief Entity-wise local assembly of a matrix from local matrices
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
 */

template <typename TMPMATRIX, class ENTITY_MATRIX_PROVIDER>
TMPMATRIX AssembleMatrixLocally(
    dim_t codim, const DofHandler &dof_handler,
    ENTITY_MATRIX_PROVIDER &entity_matrix_provider) {
  return AssembleMatrixLocally<TMPMATRIX, ENTITY_MATRIX_PROVIDER>(
      codim, dof_handler, dof_handler, entity_matrix_provider);
}

/**
 * @brief entity-local assembly of (right-hand-side) vectors from element
 * vectors
 *
 * @tparam VECTOR a generic vector type with component access through []
 * @tparam ENTITY_VECTOR_PROVIDER type for objects computing entity-local
 * vectors
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
 * - ENTITY_VECTOR_PROVIDER must
 * + offer an `Eval()` method that returns an element vector.
 * + supply an `isActive()` method for selecting cells to be taken into account
 * in assembly
 * + provide a type `ElemVec` suitable for holding an element vector
 *
 * @note Contributions of element vectors are added to the entries of the
 *       `resultvector` argument. This means that `resultvector` has to be
 *       initialized before calling this function!
 */
template <typename VECTOR, class ENTITY_VECTOR_PROVIDER>
void AssembleVectorLocally(dim_t codim, const DofHandler &dof_handler,
                           ENTITY_VECTOR_PROVIDER &entity_vector_provider,
                           VECTOR &resultvector) {
  // Pointer to underlying mesh
  auto mesh = dof_handler.Mesh();

  // Central assembly loop over entities of the co-dimension specified via
  // the template argument CODIM
  for (const lf::mesh::Entity &entity : mesh->Entities(codim)) {
    // Some cells may be skipped
    if (entity_vector_provider.isActive(entity)) {
      // Length of element vector
      const size_type veclen = dof_handler.NoLocalDofs(entity);
      // global dof indices for contribution of the entity
      nonstd::span<const gdof_idx_t> dof_idx(
          dof_handler.GlobalDofIndices(entity));
      // Request local vector from entity_vector_provider object. In the case
      // CODIM = 0, when `entity` is a cell, this is the element vector
      const auto elem_vec{entity_vector_provider.Eval(entity)};
      LF_ASSERT_MSG(elem_vec.size() >= veclen,
                    "length mismatch " << elem_vec.size() << " <-> " << veclen
                                       << ", entity " << mesh->Index(entity));
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
  VECTOR resultvector{dof_handler.NoDofs()};
  // Initialize to zero: assembly of new vector
  resultvector.setZero();
  // Perform actual assembly
  AssembleVectorLocally<VECTOR, ENTITY_VECTOR_PROVIDER>(
      codim, dof_handler, entity_vector_provider, resultvector);
  return resultvector;
}  // end AssembleVectorLocally

}  // namespace lf::assemble

#endif

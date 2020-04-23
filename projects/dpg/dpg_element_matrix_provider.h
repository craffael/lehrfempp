#ifndef PROJECTS_DPG_DPG_ELEMENT_MATRIX_PROVIDER
#define PROJECTS_DPG_DPG_ELEMENT_MATRIX_PROVIDER

/**
 * @file
 * @brief Class providing element matrices of DPG methods
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/mesh/mesh.h>

#include "product_element_matrix_provider.h"

#include "dpg.h"

namespace projects::dpg {

/**
 * @brief Class to evaluate element matrices for DPG methods.
 * @tparam SCALAR type for the entries of the element matrix. Usually 'double'.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * This class evaluates the DPG element matrix associated to a practical DPG
 * method
 *
 *
 * In an abstract setting to evaluate the DPG  element matrix the following
 * ingredients are needed
 *
 * - A diescrete Test space \f$ U_h \subset  U \f$ with
 * - An enriched discrete test space \f$ V^r \subset V \f$
 * - A bilinear form  \f$ b : U \times  V \rightarrow \mathbb{R} \f$
 * - A localizable inner product \f$ (\cdot , \cdot )_V : V \times V \rightarrow
 * \mathbb{R} \f$
 *
 * The practical DPG method can be implemented by computing a DPG element matrix
 * . If \f$ \{ b^1_u \dots b^{Q_U}_u \} \f$ denote the local shape functions of
 * \f$ U_h \f$ and \f$ \{ b^2_v \dots b^{Q_V}_v \} \f$ denote the local shape
 * functiosn of \f$ V^R \f$ on an elment \f$ K \f$, the computation of the DPG
 * element matrix boils down to
 *
 * - Evaluate the extended element stiffness Matrix \f$ B_{i,j}   =
 * b(b^j_u,b^i_v) \f$
 * - Evaluate  the local Gramian matrix \f$ G_{i,j} = (b^i_v,b^j_v)_V \f$
 * - Evaluate the DPG elment matrix \f$ A^K  = B^T G^{-1} B \f$
 */
template <typename SCALAR>
class DpgElementMatrixProvider {
 public:
  /** @brief internal type for element matrices*/
  using elem_mat_t = typename ProductElementMatrixProvider<SCALAR>::elem_mat_t;
  /** Return type of the @ref Eval() method */
  using ElemMat = typename ProductElementMatrixProvider<SCALAR>::ElemMat;

  /** @brief standard constructors */
  DpgElementMatrixProvider(const DpgElementMatrixProvider&) = delete;
  DpgElementMatrixProvider(DpgElementMatrixProvider&&) noexcept = default;
  DpgElementMatrixProvider& operator=(const DpgElementMatrixProvider&) = delete;
  DpgElementMatrixProvider& operator=(DpgElementMatrixProvider&&) = delete;

  /**
   * @brief main constructor
   * @param extendedStiffnessMatrixProvider
   * evaluates the extended element stiffness matrix \f$ B \f$
   * @param gramianProvider  evaluates the
   * local Gramian matrix \f$ G \f$
   */
  DpgElementMatrixProvider(
      std::shared_ptr<ProductElementMatrixProvider<SCALAR>>
          extendedStiffnessMatrixProvider,
      std::shared_ptr<ProductElementMatrixProvider<SCALAR>> gramianProvider)
      : extendedStiffnessMatrixProvider_(
            std::move(extendedStiffnessMatrixProvider)),
        gramianProvider_(std::move(gramianProvider)) {}

  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief main routine for the computation of DPG element matrices
   * @param cell reference to the cell for which the DPG element matrix
   * should be computed
   * @return a small dense, containing the DPG element matrix.
   */
  ElemMat Eval(const lf::mesh::Entity& cell);

  virtual ~DpgElementMatrixProvider() = default;

 private:
  /** A ProductElementMatrixProvider that evaluates the extended element
   * stiffness matrix \f$B\f$ */
  std::shared_ptr<ProductElementMatrixProvider<SCALAR>>
      extendedStiffnessMatrixProvider_;
  /** A ProductElementMatrixProvider that evaluates the local Gramian \f$ G \f$
   */
  std::shared_ptr<ProductElementMatrixProvider<SCALAR>> gramianProvider_;
};

// template deduction hint
template <class PTR>
DpgElementMatrixProvider(PTR stiffness_provider, PTR gramianProvider)
    ->DpgElementMatrixProvider<typename PTR::element_type::SCALAR>;

// evaluation method.
template <typename SCALAR>
typename DpgElementMatrixProvider<SCALAR>::ElemMat
DpgElementMatrixProvider<SCALAR>::Eval(const lf::mesh::Entity& cell) {
  // check for nullptrs
  LF_ASSERT_MSG(extendedStiffnessMatrixProvider_ != nullptr &&
                    gramianProvider_ != nullptr,
                "nullptr error for some provider");
  // check that the method is called on an active cell
  LF_ASSERT_MSG(gramianProvider_->isActive(cell) &&
                    extendedStiffnessMatrixProvider_->isActive(cell),
                "Eval method called on inactive cell. " << cell);

  // evaluate extended stiffness matrix B and local Gramian G
  ElemMat extendedStiffnessMatrix =
      extendedStiffnessMatrixProvider_->Eval(cell);
  ElemMat gramian = gramianProvider_->Eval(cell);

  // perform some size checks.
  LF_ASSERT_MSG(gramian.rows() == gramian.cols(),
                "non quadratic gramian of size ("
                    << gramian.rows() << ", " << gramian.cols() << ") on cell "
                    << cell << "\n");
  LF_ASSERT_MSG(
      extendedStiffnessMatrix.rows() == gramian.cols(),
      "size missmatch between gramian & extended Stiffness matrix on cell "
          << cell);

  // evaluate DPG element matrix A = B^T G^-1 B
  return extendedStiffnessMatrix.transpose() *
         gramian.ldlt().solve(extendedStiffnessMatrix);
}

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_ELEMENT_MATRIX_PROVIDER

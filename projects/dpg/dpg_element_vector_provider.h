#ifndef PROJECTS_DPG_DPG_ELEMENT_VECTOR_PROVIDER
#define PROJECTS_DPG_DPG_ELEMENT_VECTOR_PROVIDER

/**
 * @file
 * @brief Class providing element vectors of DPG  methods
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/mesh/mesh.h>
#include "dpg.h"
#include "product_element_matrix_provider.h"
#include "product_element_vector_provider.h"

namespace projects::dpg {
/**
 * @brief Class to evaluate element vectors of DPG methods.
 * @tparam SCALAR type for the entries of the element matrix. Usually 'double'.
 * @note This interface class complies with the requirements for the template
 *parameter `ELEM_VEC_COMP` of the function
 lf::assemble::AssembleVectorLocally().
 *
 *
 * This class evaluates the DPG element vector associated to a practical
 * DPG method
 *
 In an abstract setting to evaluate the DPG  element vector the following
 ingredients
 * are needed
 *
 * - A diescrete Test space \f$ U_h \subset  U \f$ with
 * - An enriched discrete test space \f$ V^r \subset V \f$
 * - A bilinear form  \f$ b : U \times  V \rightarrow \mathbb{R} \f$
 * - A localizable inner product \f$ (\cdot , \cdot )_V : V \times V \rightarrow
 \mathbb{R} \f$
 * - A linear form \f$ \ell: V \rightarrow  \mathbb{R} \f$
 *
 * The practical DPG method can be implemented by computing a DPG element
 vector.
 * If \f$ \{ b^1_u \dots b^{Q_U}_u \} \f$ denote the local shape functions of
 \f$ U_h \f$ and
 * \f$ \{ b^2_v \dots b^{Q_V}_v \} \f$ denote the local shape functiosn of \f$
 V^R \f$ on an elment \f$ K \f$, the
 * computation of the DPG element vector boils down to
 *
 * - Evaluate the extended element stiffness Matrix \f$ B_{i,j}   =
 b(b^j_u,b^i_v) \f$
 * - Evaluate the local Gramian matrix \f$ G_{i,j} = (b^i_v,b^j_v)_V \f$
 * - Evaluate the extended element load vector \f$ l_{i} = \ell(b^i_v) \f$
 * - Evaluate the DPG element vector \f$ \phi^K = B^T G^{-1} l \f$
 */
template <typename SCALAR>
class DpgElementVectorProvider {
 public:
  /** @brief internal type for element vectors */
  using elem_vec_t = typename ProductElementVectorProvider<SCALAR>::elem_vec_t;
  /** @brief Return type of the @ref Eval() method */
  using ElemVec = typename ProductElementVectorProvider<SCALAR>::ElemVec;
  /** @brief return type  of the Eval() method of the
   * ProductElementMatrixProvider class */
  using ElemMat = typename ProductElementMatrixProvider<SCALAR>::ElemMat;

  /** @brief standard constructors */
  DpgElementVectorProvider(const DpgElementVectorProvider&) = delete;
  DpgElementVectorProvider(DpgElementVectorProvider&&) noexcept = default;
  DpgElementVectorProvider& operator=(const DpgElementVectorProvider&) = delete;
  DpgElementVectorProvider& operator=(DpgElementVectorProvider&&) = delete;

  /**
   * @brief main constructor
   * @param extendedLoadVectorProvider
   * evaluates the extended element load vector \f$ l \f$
   * @param extendedStiffnessMatrixProvider
   * evaluates the extended element stiffness matrix \f$B \f$
   * @param gramianProvider evaluates the
   * local Gramian \f$ G \f$
   */
  DpgElementVectorProvider(
      std::shared_ptr<ProductElementVectorProvider<SCALAR>>
          extendedLoadVectorProvider,
      std::shared_ptr<ProductElementMatrixProvider<SCALAR>>
          extendedStiffnessMatrixProvider,
      std::shared_ptr<ProductElementMatrixProvider<SCALAR>> gramianProvider)
      : extendedLoadVectorProvider_(std::move(extendedLoadVectorProvider)),
        extendedStiffnessMatrixProvider_(
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
   * @brief main routine for the computation of DPG element vectors
   * @param cell refernce to the cell for which the DPG element vector
   * should be computed
   * @return a small vector, containing the DPG element vector.
   */
  ElemVec Eval(const lf::mesh::Entity& cell);

  virtual ~DpgElementVectorProvider() = default;

 private:
  /** A ProductElementVectorProvider, that evaluates the extended load vector l
   */
  std::shared_ptr<ProductElementVectorProvider<SCALAR>>
      extendedLoadVectorProvider_;
  /** A ProductElementMatrixProvider, that evaluates the extended element
   * stiffness matrix B */
  std::shared_ptr<ProductElementMatrixProvider<SCALAR>>
      extendedStiffnessMatrixProvider_;
  /** A ProductElementMatrixProvider, that evaluates the local Gramian G */
  std::shared_ptr<ProductElementMatrixProvider<SCALAR>> gramianProvider_;
};

// template deduction hint:
template <class PTR1, class PTR2>
DpgElementVectorProvider(PTR1 eLoadVectorProvider, PTR2 eStiffnessProvider,
                         PTR2 gramianProvider)
    ->DpgElementVectorProvider<typename PTR1::element_type::SCALAR>;

// evaluation method
template <typename SCALAR>
typename DpgElementVectorProvider<SCALAR>::ElemVec
DpgElementVectorProvider<SCALAR>::Eval(const lf::mesh::Entity& cell) {
  // check for nullptrs
  LF_ASSERT_MSG(extendedLoadVectorProvider_ != nullptr &&
                    extendedStiffnessMatrixProvider_ != nullptr &&
                    gramianProvider_ != nullptr,
                "nullptr error for some provider");
  LF_ASSERT_MSG(extendedLoadVectorProvider_->isActive(cell) &&
                    extendedStiffnessMatrixProvider_->isActive(cell) &&
                    gramianProvider_->isActive(cell),
                "Eval method called on inactive cell");

  // evaluate extended element vector l, extended stiffness matrix B and local
  // Gramian G.
  ElemVec extendedLoadVector = extendedLoadVectorProvider_->Eval(cell);
  ElemMat extendedStiffnessMatrix =
      extendedStiffnessMatrixProvider_->Eval(cell);
  ElemMat gramian = gramianProvider_->Eval(cell);

  // perform some size checks.
  LF_ASSERT_MSG(gramian.rows() == gramian.cols(),
                "non quadratic gramian of size ("
                    << gramian.rows() << ", " << gramian.cols() << ") on cell "
                    << cell << "\n");

  // evaluate element vector B^T G^-1 l
  return extendedStiffnessMatrix.transpose() *
         gramian.ldlt().solve(extendedLoadVector);
}
}  // namespace projects::dpg
#endif  // PROJECTS_DPG_DPG_ELEMENt_VECTOR_PROVIDER

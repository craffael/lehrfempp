#ifndef PROJECTS_DPG_PRODUCT_ELEMENT_MATRIX_PROVIDER
#define PROJECTS_DPG_PRODUCT_ELEMENT_MATRIX_PROVIDER

/**
 * @file
 * @brief Class providing element matrices associated
 * with bilinear forms between cartesian/product spaces.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <iostream>
#include <vector>

#include <lf/base/base.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"
#include "product_fe_space.h"
#include "sub_element_matrix_provider.h"

namespace projects::dpg {

/**
 * @brief Class providing element matrices associated
 * with bilinear forms between cartesian/product spaces.
 * @tparam SCALAR type for the entries of the element matrices. Ususally
 * 'double'
 * @note This interface class complies with the type requirements for the
 * template argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * This class provides element matrices associated to bilinear  forms \f$ b: U
 * \times V \rightarrow \mathbb{R} \f$  between two product spaces (c.f.
 * ProductUniformFEDofHandler, ProductUniformFESpace)
 *
 * \f[ U = U_0 \times U_1 \times \dots \times U_{n-1} \f]
 * \f[ V = V_0 \times V_1 \times \dots \times V_{m-1} \f]
 *
 * that have the following structure
 *
 * \f[ b((u_1, \dots, u_{n-1}),(v_1, \dots v_{m-1})) = \sum_{k}
 * b_k(u_{i_k},v_{j_k}) \f]
 *
 * with
 * \f[ b_k : U_{i_k} \times V_{j_k} \rightarrow \mathbb{R} \f]
 *
 * The evaluation of element matrices for such "building block" bilinear forms
 * \f$ b_k \f$ is preformed  by means of the SubElementMatrixProvider interface.
 * This class takes a vector of such
 * SubElementMatrixProviders in its constructors and builds the element matrix
 * for the bilinear form \f$ b\f$ by stacking together the sub matrices (using
 * the rules of local dof ordering specified in ProductUniformFEDofHandler).
 *
 *
 */
template <typename SCALAR>
class ProductElementMatrixProvider {
 public:
  /** @brief internal type for element matrices*/
  using elem_mat_t = typename SubElementMatrixProvider<SCALAR>::elem_mat_t;
  /** Return type of the @ref Eval() method */
  using ElemMat = typename SubElementMatrixProvider<SCALAR>::ElemMat;

  /** @brief standard constructors*/
  ProductElementMatrixProvider(const ProductElementMatrixProvider&) = delete;
  ProductElementMatrixProvider(ProductElementMatrixProvider&&) noexcept =
      default;
  ProductElementMatrixProvider& operator=(const ProductElementMatrixProvider&) =
      delete;
  ProductElementMatrixProvider& operator=(ProductElementMatrixProvider&&) =
      delete;

  /**
   * @brief main constructor
   * @param fe_space_trial the product trial finite element space  \f$ U \f$ of
   * the bilinar form
   * @param fe_space_test the product test finite element space  \f$ V \f$ of
   * the bilinear form
   * @param subproviders a vector of SubElementMatrixProviders, that evaluate
   * the element matrices associated to the  "building block" bilinear forms \f$
   * b_K \f$
   */
  ProductElementMatrixProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      std::vector<std::shared_ptr<SubElementMatrixProvider<SCALAR>>>
          subproviders)
      : fe_space_trial_(std::move(fe_space_trial)),
        fe_space_test_(std::move(fe_space_test)),
        subproviders_(std::move(subproviders)) {}

  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief main routine for the computation of  element matrices
   * @param cell refernce to the cell for which the  element matrix
   * should be computed
   * @return a small dense, containing the element matrix.
   */
  ElemMat Eval(const lf::mesh::Entity& cell);

  /** @brief virtual destructor */
  virtual ~ProductElementMatrixProvider() = default;

 private:
  /** shared pointer to the trial (product) finite element space */
  std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial_;
  /** shared pointer to the test (product) finite element space */
  std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test_;
  /** vector of pointers to providers that provide "building block" element
   * matrices */
  std::vector<std::shared_ptr<SubElementMatrixProvider<SCALAR>>> subproviders_;
};

// template deduction hint
template <typename ptr, typename vector>
ProductElementMatrixProvider(ptr fe_space_tiral, ptr fe_space_test,
                             vector subproviders)
    ->ProductElementMatrixProvider<typename ptr::element_type::SCALAR>;

// evaluation method
template <typename SCALAR>
typename ProductElementMatrixProvider<SCALAR>::ElemMat
ProductElementMatrixProvider<SCALAR>::Eval(const lf::mesh::Entity& cell) {
  // retrive dofhandlers (for local assembly with the LocalStartIndex method
  LF_ASSERT_MSG(fe_space_trial_ != nullptr && fe_space_test_ != nullptr,
                "missing description of trial or test space");
  const ProductUniformFEDofHandler& dof_h_trial{fe_space_trial_->LocGlobMap()};
  const ProductUniformFEDofHandler& dof_h_test{fe_space_test_->LocGlobMap()};

  // initialize the element matrix:
  elem_mat_t mat(dof_h_test.NoLocalDofs(cell), dof_h_trial.NoLocalDofs(cell));
  mat.setZero();

  // build the element matrix from the "building block" sub element matrizes
  for (const auto& provider : subproviders_) {
    LF_ASSERT_MSG(provider != nullptr,
                  "invalid submatrix provider, is nullptr");
    if (provider->isActive(cell)) {
      // extract the involved components
      size_type trial_component = provider->TrialComponent();
      size_type test_component = provider->TestComponent();
      // evaluate the submatrix
      elem_mat_t sub_mat = provider->Eval(cell);
      // add it to the correct place in the element matrix.
      mat.block(dof_h_test.LocalStartIndex(cell, test_component),
                dof_h_trial.LocalStartIndex(cell, trial_component),
                sub_mat.rows(), sub_mat.cols()) += sub_mat;
    }
  }
  return mat;
}

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_PRODUCT_ELEMENT_MATRIX_PROVIDER

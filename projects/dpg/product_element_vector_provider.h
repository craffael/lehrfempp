#ifndef PROJECTS_DPG_PRODUCT_ELEMENT_VECTOR_PROVIDER
#define PROJECTS_DPG_PRODUCT_ELEMENT_VECTOR_PROVIDER

/**
 * @file
 * @brief Class providing element vectors associated
 * with a linear form on a cartesian/product spaces.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <iostream>
#include <vector>

#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"
#include "product_fe_space.h"
#include "sub_element_vector_provider.h"

namespace projects::dpg {

/** type for vector length/matrix sizes */
using size_type = lf::uscalfe::size_type;

/**
 * @brief Class providing element vectors associated
 * with linear forms on cartesian/product spaces.
 * @tparam SCALAR type for the entries of the element matrices. Ususally
 * 'double'
 * @note This interface class complies with the requirements for the template
 * parameter `ELEM_VEC_COMP` of the function AssembleVectorLocally().
 *
 * This class provides element vectors associated to linar forms \f$ l: V
 * \rightarrow \mathbb R \f$  on a product space (c.f.
 * ProductUniformFEDofHandler, ProductUniformFESpace)
 *
 * \f[ V = V_0 \times V_1 \times \dots \times V_{m-1} \f]
 *
 * that has the following structure
 *
 * \f[ l((v_1, \dots v_{m-1})) = \sum_{k} l_k(v_{j_k}) \f]
 *
 * with
 * \f[ l_k :  V_{j_k} \rightarrow \mathbb{R} \f]
 *
 * The evaluation of element vectors for such "building block" linear forms \f$
 * l_k \f$ is performed by means of the SubElementVectorProvider interface.
 *
 * This class takes a vector of such SubElementVectorProviders in its
 * constructor and builds the element vector for the bilinear form \f$ l \f$ by
 * stacking together the sub vectors (using the rules of local dof ordering
 * specified in ProductUniformFEDofHandler).
 *
 */
template <typename SCALAR>
class ProductElementVectorProvider {
 public:
  /** @brief internal type for element vectors*/
  using elem_vec_t = typename SubElementVectorProvider<SCALAR>::elem_vec_t;
  /** Return type of the @ref Eval() method */
  using ElemVec = typename SubElementVectorProvider<SCALAR>::ElemVec;

  /** @brief Standard constructors*/
  ProductElementVectorProvider(const ProductElementVectorProvider&) = delete;
  ProductElementVectorProvider(ProductElementVectorProvider&&) noexcept =
      default;
  ProductElementVectorProvider& operator=(const ProductElementVectorProvider&) =
      delete;
  ProductElementVectorProvider& operator=(ProductElementVectorProvider&&) =
      delete;

  /**
   * @brief main constructor
   * @param fe_space_test the product test finite element space  \f$ V \f$ of
   * the linear form
   * @param subproviders a vector of SubElementVectorProviders, that evaluate
   *  element vectors associated to the "building block" liinear forms \f$ l_k
   * \f$
   */
  ProductElementVectorProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      std::vector<std::shared_ptr<SubElementVectorProvider<SCALAR>>>
          subproviders)
      : fe_space_test_(std::move(fe_space_test)),
        subproviders_(std::move(subproviders)) {}

  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief main routine for the computation of  element vectors
   * @param cell refernce to the cell for which the  element vector
   * should be computed
   * @return a small vector, containing the element vector.
   */
  ElemVec Eval(const lf::mesh::Entity& cell);

  /** @brief virtual destructor*/
  virtual ~ProductElementVectorProvider() = default;

 private:
  /** shared pointer to the test (product) finite element space */
  std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test_;
  /** vector of pointers to providers that provide "building block" element
   * vectors */
  std::vector<std::shared_ptr<SubElementVectorProvider<SCALAR>>> subproviders_;
};

// template deduction hint
template <typename PTR, typename vector>
ProductElementVectorProvider(PTR fe_space_test, vector subproviders)
    ->ProductElementVectorProvider<typename PTR::element_type::SCALAR>;

// evaluation method
template <typename SCALAR>
typename ProductElementVectorProvider<SCALAR>::ElemVec
ProductElementVectorProvider<SCALAR>::Eval(const lf::mesh::Entity& cell) {
  // retrive dofhandler(for local assembly with LocalStartIndex method
  LF_ASSERT_MSG(fe_space_test_ != nullptr,
                "missing description of the trial space");
  const ProductUniformFEDofHandler& dof_h{fe_space_test_->LocGlobMap()};

  // initialize the element vector
  elem_vec_t vec(dof_h.NoLocalDofs(cell));
  vec.setZero();

  // build the building block from the "building block" subproviders
  for (const auto& provider : subproviders_) {
    LF_ASSERT_MSG(provider != nullptr,
                  "invalid subvector provider, is nullptr");
    if (provider->isActive(cell)) {
      // extract involved component
      size_type test_component = provider->TestComponent();
      // evaluate the subvector
      elem_vec_t sub_vec = provider->Eval(cell);
      // add it to the corresponding place in the element vector.
      vec.segment(dof_h.LocalStartIndex(cell, test_component),
                  sub_vec.size()) += sub_vec;
    }
  }
  return vec;
}

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_PRODUCT_ELEMENT_VECTOR_PROVIDER.

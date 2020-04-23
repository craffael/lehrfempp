#ifndef PROJECTS_DPG_PRODUCT_ELEMENT_VECTOR_PROVIDER_BUILDER
#define PROJECTS_DPG_PRODUCT_ELEMENT_VECTOR_PROVIDER_BUILDER

#include "dpg.h"
#include "loc_comp_dpg.h"
#include "product_element_vector_provider.h"

namespace projects::dpg {

/**
 * @brief Builder class to build a ProductElementVectorProvider
 *
 * @tparam SCLAR type of entries of the element vectors. Fieldtype such as
 * double.
 *
 *  This class can be used to construct a ProductElementVectorProvider on a
 * product space \f[ V = V_0 \times V_1 \times \dots \times V_{m-1} \f]
 *
 * passed in the constructor. Methods are provided to add several linear forms
 * \f[ l_k :  V_{j_k} \rightarrow \mathbb{R} \f]
 * on a specified component of the space.
 *
 * the constructed ProductElementVectorProvider evalautes element vectors of the
 * linear form \f$ l: V \rightarrow \mathbb{R} \f$ given by
 *
 *  \f[ l((v_1, \dots v_{m-1})) = \sum_{k} l_k(v_{j_k}) \f]
 */
template <typename SCALAR>
class ProductElementVectorProviderBuilder {
 public:
  /** @brief standard constructor */
  ProductElementVectorProviderBuilder(
      const ProductElementVectorProviderBuilder&) = delete;
  ProductElementVectorProviderBuilder(
      ProductElementVectorProviderBuilder&&) noexcept = delete;
  ProductElementVectorProviderBuilder& operator=(
      const ProductElementVectorProviderBuilder&) = delete;
  ProductElementVectorProviderBuilder& operator=(
      ProductElementVectorProviderBuilder&&) = delete;

  /**
   * @brief main contructor, construct a new builder
   * @param fe_space_test collection of specifications about the  (product)
   * space \f$ V \f$
   */
  explicit ProductElementVectorProviderBuilder(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test)
      : fe_space_test_(std::move(fe_space_test)) {}

  /**
   * @brief Adds a linear form \f$ l_k \f$ to \f$ l \f$ representing a simple
   * load linear form.
   * @param test_component index of the  component  \f$ v \f$ of \f$ l_k \f$
   * @param f mesh function for the scalar valued source function \f$ f \f$
   *
   * @tparam FUNCTOR see the type requirements of the template parameter FUNCTOR
   * of the LoadElementVectorProvider class.
   *
   * The (local)  added linear form \f$ l_k\f$   is
   @f[
      v \mapsto \int_K f(\mathbf{x})\,v\,\mathrm{d}\mathbf{x}\;,
 * @f]
   *
   * For further information about the added  added linear form \f$ l_k \f$ see
   * the documentation of
   * LoadElementVectorProvider
   */
  template <typename FUNCTOR>
  ProductElementVectorProviderBuilder& AddLoadElementVectorProvider(
      size_type test_component, FUNCTOR f);

  /**
   * @brief Build the ProductElementVectorProvider based on the provided linear
   * forms
   */
  std::shared_ptr<ProductElementVectorProvider<SCALAR>> Build();

  /** @brief default destructor */
  ~ProductElementVectorProviderBuilder() = default;

 private:
  /** @brief collection of specifications for the test space */
  std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test_;
  /** @brief vector of Element Vector Providers representing already added
   * linear forms */
  std::vector<std::shared_ptr<SubElementVectorProvider<SCALAR>>> subproviders_;
};

template <typename SCALAR>
template <typename FUNCTOR>
ProductElementVectorProviderBuilder<SCALAR>&
ProductElementVectorProviderBuilder<SCALAR>::AddLoadElementVectorProvider(
    size_type test_component, FUNCTOR f) {
  subproviders_.push_back(
      std::make_shared<LoadElementVectorProvider<SCALAR, FUNCTOR>>(
          fe_space_test_, test_component, f));
  return *this;
}

template <typename SCALAR>
std::shared_ptr<ProductElementVectorProvider<SCALAR>>
ProductElementVectorProviderBuilder<SCALAR>::Build() {
  auto provider = std::make_shared<ProductElementVectorProvider<SCALAR>>(
      fe_space_test_, subproviders_);

  // clear supplied information
  subproviders_ = {};
  return provider;
}

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_PRODUCT_ELEMENT_VECTOR_PROVIDER_BUILDER

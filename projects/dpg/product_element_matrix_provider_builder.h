#ifndef PROJECTS_DPG_PRODUCT_ELEMENT_MATRIX_PROVIDER_BUILDER
#define PROJECTS_DPG_PRODUCT_ELEMENT_MATRIX_PROVIDER_BUILDER

#include "dpg.h"
#include "loc_comp_dpg.h"
#include "product_element_matrix_provider.h"

namespace projects::dpg {
/**
 * @brief Builder class to build a  ProductElementMatrixProvider
 *
 * @tparam SCALAR type of entries of the element matrics. Field type such as
 * double
 *
 * This class can be used to construct a ProductElementMatrixProvider between
 * two product spaces
 *
 * \f[ U = U_0 \times U_1 \times \dots \times U_{n-1} \f]
 * \f[ V = V_0 \times V_1 \times \dots \times V_{m-1} \f]
 *
 * passed in the constructor. Methods are provided to add several bilinear forms
 * \f[ b_k : U_{i_k} \times V_{j_k} \rightarrow \mathbb{R} \f]
 * between two specified components of the spaces to the bilinear form.
 *
 * The constructed  ProductElementMatrixProvider evaluates element matrices of
 * the bilinear form \f$ b: U \times V \rightarrow \mathbb{R} \f$ given by
 *
 * \f[ b((u_1, \dots, u_{n-1}),(v_1, \dots v_{m-1})) = \sum_{k}
 * b_k(u_{i_k},v_{j_k}) \f]
 */
template <typename SCALAR>
class ProductElementMatrixProviderBuilder {
 public:
  /** @brief standard constructor */
  ProductElementMatrixProviderBuilder(
      const ProductElementMatrixProviderBuilder&) = delete;
  ProductElementMatrixProviderBuilder(
      ProductElementMatrixProviderBuilder&&) noexcept = delete;
  ProductElementMatrixProviderBuilder& operator=(
      const ProductElementMatrixProviderBuilder&) = delete;
  ProductElementMatrixProviderBuilder& operator=(
      ProductElementMatrixProviderBuilder&&) = delete;

  /**
   * @brief main constructor, constructs a new builder
   * @param fe_space_trial collection of specifications about the  (product)
   * trial space \f$ U \f$
   * @param fe_space_test collection of specifications about the  (product) test
   * space \f$ V \f$
   */
  ProductElementMatrixProviderBuilder(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test)
      : fe_space_trial_(std::move(fe_space_trial)),
        fe_space_test_(std::move(fe_space_test)) {}

  /**
   * @brief Adds a bilinear form \f$ b_k \f$ to \f$ b \f$ representing a
 diffusion term.
   * @param trial_component index of the trial component \f$ u \f$ of \f$ b_k
 \f$
   * @param test_component index of the test component  \f$ v \f$ of \f$ b_k \f$
   * @param alpha mesh function of the matrix or scalar valued diffusion
   coefficient \f$ \alpha \f$
   *
   * @tparam DIFF_COEFF see the type requirements of the template parameter
   DIFF_COEFF
   * in the DiffusionElementMatrixProvider class.
   *
   * The (local)  added bilinear form \f$ b_k\f$   is
   *
 @f[
    (u,v) \mapsto\int\limits_{K}\mathbf{grad}\,u
          \cdot \boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,v
 \mathrm{d}\mathbf{x}
 \;,
 * @f]
   *
   * For further information about the added bilinear form \f$ b_k \f$
   * see documentation of the class DiffusionElementMatrixProvider
   *
   */
  template <typename DIFF_COEFF>
  ProductElementMatrixProviderBuilder& AddDiffusionElementMatrixProvider(size_type trial_component,
                                              size_type test_component,
                                              DIFF_COEFF alpha);

  /**
   * @brief Adds a bilinear form \f$ b_k \f$ to \f$ b \f$ representing a
 convection-like term.
   * @param trial_component index of the trial component \f$ u \f$ of \f$ b_k
 \f$
   * @param test_component index of the test component  \f$ v \f$ of \f$ b_k \f$
   * @param beta_1 mesh function for the "first" vector valued convection
   coefficient \f$ \beta_1 \f$
   * @param beta_2 mesh function for the "second" vector valued convection
   coefficient \f$ \beta_2 \f$
   *
   * @tparam CONVECTION_COEFF_1 and CONVECTION_COEFF_2 see the type requirements
   of the template paramters
   * CONVECTION_COEFF_1 and CONVECTION_COEFF_2 in the
   ConvectionElementMatrixProvider class.
   *
   *  The (local)  added bilinear form \f$ b_k\f$   is
   *
   * @f[
    (u,v) \mapsto\int\limits_{K}
        \mathbf{grad}\,v \cdot \boldsymbol{\beta_1}(\mathbf{x})u +
        \mathbf{grad}\,u \cdot \boldsymbol{\beta_2}(\mathbf{x})v
 \mathrm{d}\mathbf{x}
 \;,
 * @f]
   *
   * For further information about the added bilinear form \f$ b_k \f$  see
   * the documentation of ConvectionElementMatrixProvider
   */
  template <typename CONVECTION_COEFF_1, typename CONVECTION_COEFF_2>
  ProductElementMatrixProviderBuilder& AddConvectionElementMatrixProvider(size_type trial_component,
                                               size_type test_component,
                                               CONVECTION_COEFF_1 beta_1,
                                               CONVECTION_COEFF_2 beta_2);

  /**
   * @brief Adds a bilinear form \f$ b_k \f$ to \f$ b \f$ representing a
 reaction term.
   * @param trial_component index of the trial component \f$ u \f$ of \f$ b_k
 \f$
   * @param test_component index of the test component  \f$ v \f$ of \f$ b_k \f$

   * @param gamma mesh function for the scalar valued reaction coefficient \f$
 \gamma \f$.
   *
   * @tparam REACTION_COEFF see the type requirements for the template parameter
   REACTION_COEFF
   * in the ReactionElementMatrixProvider class.
   *
   *  The (local)  added bilinear form \f$ b_k\f$   is
 * @f[
    (u,v) \mapsto\int\limits_{K} \gamma(\mathbf{x})u\,v\,\mathrm{d}\mathbf{x}
 \;,
 * @f]
   *
   * For further information about the added bilinear form \f$ b_k \f$ see
   * the documentation of ReactionElementMatrixProvider.
   */
  template <typename REACTION_COEFF>
  ProductElementMatrixProviderBuilder& AddReactionElementMatrixProvider(size_type trial_component,
                                             size_type test_component,
                                             REACTION_COEFF gamma);

  /**
   * @brief Adds a bilinear form \f$ b_k \f$ to \f$ b \f$ representing a term
 that involves
   * a fulx component in the trial space.

   * @param trial_component index of the trial component \f$ \hat{q}_n \f$ of
 \f$ b_k \f$.
   * This component should represent a flux.
   @param test_component index of the test component  \f$ v \f$ of \f$ b_k \f$

   * @param alpha mesh function for the scalar valued coefficient \f$ \alpha \f$
   *
   * @tparam DIFF_COEFF see the type requeirements for the template parameter
   DIFF_COEFF
   * in the FluxElementMatrixProvider class.
   *
   * The (local)  added bilinear form \f$ b_k\f$   is
   *
  * @f[
    (\hat{q}_n,v) \mapsto\int\limits_{\partial K} \hat{q}_n
 \mathrm{sgn}_K(\mathbf{x}) \alpha(\mathbf{x}) v \mathrm{d}\mathbf{x}
 \;,
  * @f]
   *
   * For further information about the added bilinear form \f$ b_k \f$
   *  see the documentation of FluxElementMatrixProvider.
   */
  template <typename DIFF_COEFF>
  ProductElementMatrixProviderBuilder& AddFluxElementMatrixProvider(size_type trial_component,
                                         size_type test_component,
                                         DIFF_COEFF alpha);


  /**
   * @brief Adds a bilinear form \f$ b_k \f$ to \f$ b \f$ representing a term
 that involves
   * a trace component in the trial space.

   * @param trial_component index of the trial component \f$ \hat{u} \f$ of \f$
 b_k \f$.
   * This component should represent a trace.
   @param test_component index of the test component  \f$ v \f$ of \f$ b_k \f$
   * @param beta mesh function for the vector valued coefficient \f$ \beta \f$
   *
   * @tparam COEFF see the type requirements for the template parameter COEFF in
   * the TraceElementMAtrixProvider class.
   *
   * The (local)  added bilinear form \f$ b_k\f$   is
 * @f[
    (\hat u,v) \mapsto\int\limits_{\partial K} \hat u \mathbf{n}_K \cdot
 \boldsymbol{\beta}(\mathbf{x}) v\mathrm{d}\mathbf{x}
 \;,
 * @f]
   * For further information about the added bilinear form \f$ b_k \f$
   *  see the documentation of TraceElementMatrixProvider.
   */
  template <typename COEFF>
  ProductElementMatrixProviderBuilder& AddTraceElementMatrixProvider(size_type trial_component,
                                          size_type test_component, COEFF beta);

  /**
   * @brief Build the ProdcutElementMatrixProvider for the bilinear form \f$
   * b\f$ based on the provided bilinear forms \f$ b_k \f$
   */
  std::shared_ptr<ProductElementMatrixProvider<SCALAR>> Build();

  /** @brief virtual destructor */
  ~ProductElementMatrixProviderBuilder() = default;

 private:
  /** @brief  collection of specifications for the trial space */
  std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial_;
  /** @brief collection of specifications for the test space */
  std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test_;
  /** @brief vector of Element Matrix Providers representing already added
   * bilinear forms */
  std::vector<std::shared_ptr<SubElementMatrixProvider<SCALAR>>> subproviders_;
};

template <typename SCALAR>
template <typename DIFF_COEFF>
ProductElementMatrixProviderBuilder<SCALAR>&
ProductElementMatrixProviderBuilder<SCALAR>::AddDiffusionElementMatrixProvider(
    size_type trial_component, size_type test_component, DIFF_COEFF alpha) {
  subproviders_.push_back(
      std::make_shared<DiffusionElementMatrixProvider<SCALAR, DIFF_COEFF>>(
          fe_space_trial_, fe_space_test_, trial_component, test_component,
          alpha));

  return *this;
}

template <typename SCALAR>
template <typename CONVECTION_COEFF_1, typename CONVECTION_COEFF_2>
ProductElementMatrixProviderBuilder<SCALAR>&
ProductElementMatrixProviderBuilder<SCALAR>::AddConvectionElementMatrixProvider(
    size_type trial_component, size_type test_component,
    CONVECTION_COEFF_1 beta_1, CONVECTION_COEFF_2 beta_2) {
  subproviders_.push_back(std::make_shared<ConvectionElementMatrixProvider<
                              SCALAR, CONVECTION_COEFF_1, CONVECTION_COEFF_2>>(
      fe_space_trial_, fe_space_test_, trial_component, test_component, beta_1,
      beta_2));
  return *this;
}

template <typename SCALAR>
template <typename REACTION_COEFF>
ProductElementMatrixProviderBuilder<SCALAR>&
ProductElementMatrixProviderBuilder<SCALAR>::AddReactionElementMatrixProvider(
    size_type trial_component, size_type test_component, REACTION_COEFF gamma) {
  subproviders_.push_back(
      std::make_shared<ReactionElementMatrixProvider<SCALAR, REACTION_COEFF>>(
          fe_space_trial_, fe_space_test_, trial_component, test_component,
          gamma));
  return *this;
}

template <typename SCALAR>
template <typename DIFF_COEFF>
ProductElementMatrixProviderBuilder<SCALAR>&
ProductElementMatrixProviderBuilder<SCALAR>::AddFluxElementMatrixProvider(
    size_type trial_component, size_type test_component, DIFF_COEFF alpha) {
  subproviders_.push_back(
      std::make_shared<FluxElementMatrixProvider<SCALAR, DIFF_COEFF>>(
          fe_space_trial_, fe_space_test_, trial_component, test_component,
          alpha));
  return *this;
}

template <typename SCALAR>
template <typename COEFF>
ProductElementMatrixProviderBuilder<SCALAR>&
ProductElementMatrixProviderBuilder<SCALAR>::AddTraceElementMatrixProvider(
    size_type trial_component, size_type test_component, COEFF beta) {
  subproviders_.push_back(
      std::make_shared<TraceElementMatrixProvider<SCALAR, COEFF>>(
          fe_space_trial_, fe_space_test_, trial_component, test_component,
          beta));
  return *this;
}

template <typename SCALAR>
std::shared_ptr<ProductElementMatrixProvider<SCALAR>>
ProductElementMatrixProviderBuilder<SCALAR>::Build() {
  auto provider = std::make_shared<ProductElementMatrixProvider<SCALAR>>(
      fe_space_trial_, fe_space_test_, subproviders_);
  // clear supplied information
  subproviders_ = {};
  return provider;
}

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_PRODUCT_ELEMENT_MATRIX_PROVIDER_BUILDER

#ifndef PROJECTS_DPG_LOC_COMP_DPG
#define PROJECTS_DPG_LOC_COMP_DPG

/**
 * @file
 * @brief classes implementing the SubElementMatrixProvider and
 * SubElementVectorProvider interfaces, performing local computations.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"
#include "dpg_tools.h"
#include "product_fe_space.h"
#include "sub_element_matrix_provider.h"
#include "sub_element_vector_provider.h"

namespace projects::dpg {

/**
 * @brief data structure for passing collections of quadrature rules,
 * see lf::uscalfe for more information.
 */
using quad_rule_collection_t = lf::uscalfe::quad_rule_collection_t;

/**
 * @headerfile projects/dpg/loc_comp_dpg.h
 * @brief Class for local quadrature based computations of  sub element matrices
 corresponding
 * to  diffusion element matrices.
 *
 * @tparam SCALAR type for the entries of the element matrices. Must be a field
 *                     type such as `double` or `std::complex<double>`
 * @tparam DIFF_COEFF a \ref mesh_function "MeshFunction" that defines the
 *                    diffusion coefficient \f$ \mathbf{\alpha} \f$.
 *                    It should be either scalar- or matrix-valued.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * The element matrix corresponds to the (local) bilinear form
 * @f[
    (u,v) \mapsto\int\limits_{K}\mathbf{grad}\,u
          \cdot \boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,v
 \mathrm{d}\mathbf{x}
 \;,
 * @f]
 * with diffusion coefficient @f$\mathbf{\alpha}@f$.
 * \f$u\f$ is a component of the trial space and \f$v\f$ a component of the test
 space.
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - DIFF_COEFF must provide an evaluation operator
 * `operator (const Entity &,ref_coord_t)` that returns either a scalar
 * or a matrix type that is compatible with Eigen's matrices. Usually it will
 * be an Eigen::Matrix either of variable or fixed size.
 *
 */
template <typename SCALAR, typename DIFF_COEFF>
class DiffusionElementMatrixProvider : public SubElementMatrixProvider<SCALAR> {
  static_assert(lf::uscalfe::isMeshFunction<DIFF_COEFF>);

 public:
  /** @brief inherited types for element matrices */
  using elem_mat_t = typename SubElementMatrixProvider<SCALAR>::elem_mat_t;
  using ElemMat = typename SubElementMatrixProvider<SCALAR>::ElemMat;

  DiffusionElementMatrixProvider(const DiffusionElementMatrixProvider&) =
      delete;
  DiffusionElementMatrixProvider(DiffusionElementMatrixProvider&&) noexcept =
      default;
  DiffusionElementMatrixProvider& operator=(
      const DiffusionElementMatrixProvider&) = delete;
  DiffusionElementMatrixProvider& operator=(DiffusionElementMatrixProvider&&) =
      delete;

  /**
   * @brief Constructor: performs cell-independent precomputations.
   * @param fe_space_trial collection of specifications for the trial fe space
   * @param fe_space_test collection of specifications for the test fe space
   * @param trial_component Index of the trial space component \f$u\f$
   * @param test_component Index of the  test space component \f$v\f$
   * @param alpha mesh function for the (matrix or scalar valued) diffusion
   * coefficient.
   *
   * This constructur uses local quadrature rules with the degree of exactness
   * chosen as the sum of the polynomial degrees of the two components \f$ u\f$
   * and \f$ v \f$
   */
  DiffusionElementMatrixProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      size_type trial_component, size_type test_component, DIFF_COEFF alpha);

  /**
   * @brief main routine for the computation of element matrices
   * @param cell reference to a (triangular or quadrilateral) cell for which the
   *        element matirx should be computed.
   * @return small dense matrix containing the element matrix
   *
   * Actual computation is based on numerical quadrature and mapping techniques.
   *
   * Throws an assertion, in case any specification is missing for the type of
   * cell.
   */
  ElemMat Eval(const lf::mesh::Entity& cell) override;

  [[nodiscard]] size_type TrialComponent() const override {
    return trial_component_;
  }

  [[nodiscard]] size_type TestComponent() const override {
    return test_component_;
  }

  ~DiffusionElementMatrixProvider() override = default;

 private:
  /** @brief functor providing the diffusion coefficient */
  DIFF_COEFF alpha_;

  /** @brief array containing the precomputed information for the trial space:
   * fe_precomp_trial[i] contains the precomputed reference finite element of
   * the trial space component u for ref_el i.*/
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_trial_;
  /** @brief array containing the precomputed information for the test space:
   * fe_precomp_trial[i] contains the precomputed reference finite element of
   * the test space component v for ref_el i.*/
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_test_;

  /** @brief index of trial component u */
  size_type trial_component_;
  /** @brief index of test component v */
  size_type test_component_;
};

// template deduction hint.
template <class PTR, class DIFF_COEFF>
DiffusionElementMatrixProvider(PTR fe_space_trial, PTR fe_space_test,
                               size_type trial_component,
                               size_type test_component, DIFF_COEFF alpha)
    ->DiffusionElementMatrixProvider<typename PTR::element_type::SCALAR,
                                     DIFF_COEFF>;

// Constructor. internal construction of quadrature rules.
template <typename SCALAR, typename DIFF_COEFF>
DiffusionElementMatrixProvider<SCALAR, DIFF_COEFF>::
    DiffusionElementMatrixProvider(
        std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
        std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
        size_type trial_component, size_type test_component, DIFF_COEFF alpha)
    : alpha_(std::move(alpha)),
      trial_component_(trial_component),
      test_component_(test_component),
      fe_precomp_trial_(),
      fe_precomp_test_() {
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    // obtain descriptions of local shape functions.
    auto fe_trial =
        fe_space_trial->ShapeFunctionLayout(ref_el, trial_component_);
    auto fe_test = fe_space_test->ShapeFunctionLayout(ref_el, test_component_);

    // check, that  shape functions for both the trial and test space component
    // for that entity are available.
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe_trial != nullptr && fe_test != nullptr) {
      // use a quadrature rule whose degree is the sum of the degrees of the fe
      // spaces.
      size_type degree = fe_trial->Degree() + fe_test->Degree();
      lf::quad::QuadRule qr = lf::quad::make_QuadRule(ref_el, degree);
      // Precompute cell-independent quantities.
      fe_precomp_trial_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(fe_trial, qr);
      fe_precomp_test_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(fe_test, qr);
    }
  }
}

template <typename SCALAR, typename DIFF_COEFF>
typename DiffusionElementMatrixProvider<SCALAR, DIFF_COEFF>::ElemMat
DiffusionElementMatrixProvider<SCALAR, DIFF_COEFF>::Eval(
    const lf::mesh::Entity& cell) {
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>& pfe_trial =
      fe_precomp_trial_[cell.RefEl().Id()];
  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>& pfe_test =
      fe_precomp_test_[cell.RefEl().Id()];

  // check initialization:
  LF_ASSERT_MSG(pfe_trial.isInitialized(),
                "No local shape function on trial space for entity type "
                    << cell.RefEl());
  LF_ASSERT_MSG(
      pfe_test.isInitialized(),
      "No local shape function on test space for entity type " << cell.RefEl());
  // check the consistency of the quadratrue rules used.
  LF_ASSERT_MSG(pfe_trial.Qr().Points() == pfe_test.Qr().Points() &&
                    pfe_trial.Qr().Weights() == pfe_test.Qr().Weights() &&
                    pfe_trial.Qr().NumPoints() == pfe_test.Qr().NumPoints(),
                "Qr missmatch");

  // query the geometry of the cell.
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry");
  LF_ASSERT_MSG(geo_ptr->DimLocal() == 2, "Only 2D implementation available");

  // phsical dimension of the cell.
  const lf::uscalfe::dim_t world_dim = geo_ptr->DimGlobal();

  // retrive information for parametric fem computation:
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe_trial.Qr().Points()));
  LF_ASSERT_MSG(determinants.size() == pfe_trial.Qr().NumPoints(),
                "Mismatch " << determinants.size() << " <-> "
                            << pfe_trial.Qr().NumPoints());

  const Eigen::MatrixXd JinvT(
      geo_ptr->JacobianInverseGramian(pfe_trial.Qr().Points()));
  LF_ASSERT_MSG(
      JinvT.cols() == 2 * pfe_trial.Qr().NumPoints() &&
          JinvT.rows() == world_dim,
      "Mismatch " << JinvT.cols() << " <-> " << 2 * pfe_trial.Qr().NumPoints());

  // evaluate alpha in the quadrature points
  auto alphaeval = alpha_(cell, pfe_trial.Qr().Points());

  // initialize the element matrix
  elem_mat_t mat(pfe_test.NumRefShapeFunctions(),
                 pfe_trial.NumRefShapeFunctions());
  mat.setZero();

  // loop over quadrature points:
  for (size_type k = 0; k < pfe_trial.Qr().NumPoints(); ++k) {
    // transformed quadrature weight.
    const double w = pfe_trial.Qr().Weights()[k] * determinants[k];

    // transformed gradients
    const auto trf_trial_grad(
        JinvT.block(0, 2 * k, world_dim, 2) *
        pfe_trial.PrecompGradientsReferenceShapeFunctions()
            .block(0, 2 * k, mat.cols(), 2)
            .transpose());
    const auto trf_test_grad(JinvT.block(0, 2 * k, world_dim, 2) *
                             pfe_test.PrecompGradientsReferenceShapeFunctions()
                                 .block(0, 2 * k, mat.rows(), 2)
                                 .transpose());

    mat += w * trf_test_grad.transpose() * alphaeval[k] * trf_trial_grad;
  }
  return mat;
}

/**
 * @headerfile projects/dpg/loc_comp_dpg.h
 * @brief Class for local quadrature based computations of sub element matrices
 corresponding to
 * reaction element matrices.
 *
 * @tparam SCALAR type for the entries of the element matrices. Must be a field
 *                     type such as `double` or `std::complex<double>`
 * @tparam REACTION_COEFF a \ref mesh_function "MeshFunction" that defines the
 *                    reaction coefficient \f$ \gamma \f$.
 *                    It should be scalar valued
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * The element matrix corresponds to the (local) bilinear form
 * @f[
    (u,v) \mapsto\int\limits_{K} \gamma(\mathbf{x})u\,v\,\mathrm{d}\mathbf{x}
 \;,
 * @f]
 *
 * with reaction coefficient \f$\gamma \f$. \f$u\f$ is a component of the trial
 space and \f$v\f$ a component of the test space.
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - REACTION_COEFF must provide an evaluation operator
 * `operator (const Entity &,ref_coord_t)` that returns a scalar type. Usually
 it will be double.
 *
 */
template <typename SCALAR, typename REACTION_COEFF>
class ReactionElementMatrixProvider : public SubElementMatrixProvider<SCALAR> {
  static_assert(lf::uscalfe::isMeshFunction<REACTION_COEFF>);

 public:
  /** @brief inherited types for element matrices */
  using elem_mat_t = typename SubElementMatrixProvider<SCALAR>::elem_mat_t;
  using ElemMat = typename SubElementMatrixProvider<SCALAR>::ElemMat;

  ReactionElementMatrixProvider(const ReactionElementMatrixProvider&) = delete;
  ReactionElementMatrixProvider(ReactionElementMatrixProvider&&) noexcept =
      default;
  ReactionElementMatrixProvider& operator=(
      const ReactionElementMatrixProvider&) = delete;
  ReactionElementMatrixProvider& operator=(ReactionElementMatrixProvider&&) =
      delete;

  /**
   * @brief Constructor: performs cell-independent precomputations.
   * @param fe_space_trial collection of specifications for the trial fe space
   * @param fe_space_test collection of  specifications for the test  fe space
   * @param trial_component Index of the trial space component \f$u\f$
   * @param test_component  Index of the  test space component \f$v\f$
   * @param gamma mesh function for the (scalar valued) reaction coefficient.
   *
   * This constructur uses local quadrature rules with the degree of exactness
   * chosen as the sum of the polynomial degrees of the two components \f$ u\f$
   * and \f$ v \f$
   */
  ReactionElementMatrixProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      size_type trial_component, size_type test_component,
      REACTION_COEFF gamma);

  /**
   * @brief main routine for the computation of element matrices
   * @param cell reference to a (triangular or quadrilateral) cell for which the
   *        element matirx should be computed.
   * @return small dense matrix containing the element matrix
   *
   * Actual computation is based on numerical quadrature and mapping techniques.
   *
   * Throws an assertion, in case any specification is missing for the type of
   * cell.
   */
  ElemMat Eval(const lf::mesh::Entity& cell) override;

  [[nodiscard]] size_type TrialComponent() const override {
    return trial_component_;
  }

  [[nodiscard]] size_type TestComponent() const override {
    return test_component_;
  }

  ~ReactionElementMatrixProvider() override = default;

 private:
  /** @brief functor providing reaction coefficient */
  REACTION_COEFF gamma_;

  /** @brief array containing the precomputed information for the trial space:
   * fe_precomp_trial[i] contains the precomputed reference finite element of
   * the trial space component u for ref_el i.*/
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_trial_;
  /** @brief array containing the precomputed information for the test space:
   * fe_precomp_trial[i] contains the precomputed reference finite element of
   * the test space component v for ref_el i.*/
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_test_;

  /** @brief index of trial component u */
  size_type trial_component_;

  /** @brief index of test component v */
  size_type test_component_;
};

// template deduction hint
template <class PTR, class REACTION_COEFF>
ReactionElementMatrixProvider(PTR fe_trial, PTR fe_test,
                              size_type trial_component,
                              size_type test_component, REACTION_COEFF gamma)
    ->ReactionElementMatrixProvider<typename PTR::element_type::SCALAR,
                                    REACTION_COEFF>;

// Constructor. internal construction of quadrature rules.
template <typename SCALAR, typename REACTION_COEFF>
ReactionElementMatrixProvider<SCALAR, REACTION_COEFF>::
    ReactionElementMatrixProvider(
        std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
        std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
        size_type trial_component, size_type test_component,
        REACTION_COEFF gamma)
    : gamma_(std::move(gamma)),
      trial_component_(trial_component),
      test_component_(test_component),
      fe_precomp_trial_(),
      fe_precomp_test_() {
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    // obtain descriptions of local shape functions.
    auto fe_trial =
        fe_space_trial->ShapeFunctionLayout(ref_el, trial_component_);
    auto fe_test = fe_space_test->ShapeFunctionLayout(ref_el, test_component_);

    // check, that  shape functions for both the trial and test space component
    // for that entity are available.
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe_trial != nullptr && fe_test != nullptr) {
      // use a quadrature rule whose degree is the sum of the degrees of the fe
      // spaces.
      size_type degree = fe_trial->Degree() + fe_test->Degree();
      lf::quad::QuadRule qr = lf::quad::make_QuadRule(ref_el, degree);
      // Precompute cell-independent quantities.
      fe_precomp_trial_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(fe_trial, qr);
      fe_precomp_test_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(fe_test, qr);
    }
  }
}

template <typename SCALAR, typename REACTION_COEFF>
typename ReactionElementMatrixProvider<SCALAR, REACTION_COEFF>::ElemMat
ReactionElementMatrixProvider<SCALAR, REACTION_COEFF>::Eval(
    const lf::mesh::Entity& cell) {
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>& pfe_trial =
      fe_precomp_trial_[cell.RefEl().Id()];
  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>& pfe_test =
      fe_precomp_test_[cell.RefEl().Id()];

  // check initialization:
  LF_ASSERT_MSG(pfe_trial.isInitialized(),
                "No local shape function on trial space for entity type "
                    << cell.RefEl());
  LF_ASSERT_MSG(
      pfe_test.isInitialized(),
      "No local shape function on test space for entity type " << cell.RefEl());
  // check the consistency of the quadratrue rules used.
  LF_ASSERT_MSG(pfe_trial.Qr().Points() == pfe_test.Qr().Points() &&
                    pfe_trial.Qr().Weights() == pfe_test.Qr().Weights() &&
                    pfe_trial.Qr().NumPoints() == pfe_test.Qr().NumPoints(),
                "Qr missmatch");

  // query the geometry of the cell.
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry");
  LF_ASSERT_MSG(geo_ptr->DimLocal() == 2, "Only 2D implementation available");

  // const lf::uscalfe::dim_t world_dim = geo_ptr->DimGlobal();

  // retrive information for parametric fem computation:
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe_trial.Qr().Points()));
  LF_ASSERT_MSG(determinants.size() == pfe_trial.Qr().NumPoints(),
                "Mismatch " << determinants.size() << " <-> "
                            << pfe_trial.Qr().NumPoints());

  // evaluate gamm in the quadrature points.
  auto gammaeval = gamma_(cell, pfe_trial.Qr().Points());

  // initialize element matrix.
  elem_mat_t mat(pfe_test.NumRefShapeFunctions(),
                 pfe_trial.NumRefShapeFunctions());
  mat.setZero();

  for (size_type k = 0; k < pfe_trial.Qr().NumPoints(); ++k) {
    // transformed quadratrue weight.
    const double w = pfe_trial.Qr().Weights()[k] * determinants[k];
    // update element matrix.
    mat += w * pfe_test.PrecompReferenceShapeFunctions().col(k) * gammaeval[k] *
           pfe_trial.PrecompReferenceShapeFunctions().col(k).transpose();
  }

  return mat;
}

/**
 * @headerfile projects/dpg/loc_comp_dpg.h
 * @brief Class for local quadrature based computations sub element matrices
 corresponding
 * to convection-like element matrices.
 *
 * @tparam SCALAR type for the entries of the element matrices. Must be a field
 *                     type such as `double` or `std::complex<double>`
 * @tparam CONVECTION_COEFF_1 a \ref mesh_function "MeshFunction" that defines
 the
 *                    "first" convection coefficient \f$ \mathbf{\beta_1} \f$.
 *                    It should be vector valued.
 * @tparam CONVECTION_COEFF_2 a \ref mesh_function "MeshFunction" that defines
 the
 *                    "second" convection coefficient \f$ \mathbf{\beta_2} \f$.
 *                    It should be vector valued.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * The element matrix corresponds to the (local) bilinear form
 * @f[
    (u,v) \mapsto\int\limits_{K}
        \mathbf{grad}\,v \cdot \boldsymbol{\beta_1}(\mathbf{x})u +
        \mathbf{grad}\,u \cdot \boldsymbol{\beta_2}(\mathbf{x})v
 \mathrm{d}\mathbf{x}
 \;,
 * @f]
 * with first "convection" coefficient @f$\mathbf{\beta_1}@f$ and second
 "convection" coefficient
 * @f$\mathbf{\beta_2}@f$.  \f$ u \f$ is a component of the trial space and \f$
 v\f$  a component of the test space.
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - CONVECTION_COEFF_1 and CONVECTION_COEFF_2 must provide an evaluation
 operator
 * `operator (const Entity &,ref_coord_t)`  that returns a vector that is
 compatible with Eigen's matrices. Usually it will
 * be an Eigen::Vector of variable or fixed size.
 *
 */
template <typename SCALAR, typename CONVECTION_COEFF_1,
          typename CONVECTION_COEFF_2>
class ConvectionElementMatrixProvider
    : public SubElementMatrixProvider<SCALAR> {
  static_assert(lf::uscalfe::isMeshFunction<CONVECTION_COEFF_1>);
  static_assert(lf::uscalfe::isMeshFunction<CONVECTION_COEFF_2>);

 public:
  /** @brief inherited types for element matrices */
  using elem_mat_t = typename SubElementMatrixProvider<SCALAR>::elem_mat_t;
  using ElemMat = typename SubElementMatrixProvider<SCALAR>::ElemMat;

  ConvectionElementMatrixProvider(const ConvectionElementMatrixProvider&) =
      delete;
  ConvectionElementMatrixProvider(ConvectionElementMatrixProvider&&) noexcept =
      default;
  ConvectionElementMatrixProvider& operator=(
      const ConvectionElementMatrixProvider&) = delete;
  ConvectionElementMatrixProvider& operator=(
      ConvectionElementMatrixProvider&&) = delete;

  /**
   * @brief Constructor: performs cell-independent precomputations.
   * @param fe_space_trial collection of specifications for the trial fe space
   * @param fe_space_test collection of specifications for the test fe space
   * @param trial_component  Index of the trial space component \f$u\f$
   * @param test_component Index of the  test space component \f$v\f$
   * @param beta_1 mesh function for the "first" vector valued convection
   * coefficient
   * @param beta_2 mesh function for the "second" vector valued convection
   * coefficient
   *
   * This constructur uses local quadrature rules with the degree of exactness
   * chosen as the sum of the polynomial degrees of the two components \f$ u\f$
   * and \f$ v \f$
   */
  ConvectionElementMatrixProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      size_type trial_component, size_type test_component,
      CONVECTION_COEFF_1 beta_1, CONVECTION_COEFF_2 beta_2);

  /**
   * @brief main routine for the computation of element matrices
   * @param cell reference to a (triangular or quadrilateral) cell for which the
   *        element matirx should be computed.
   * @return small dense matrix containing the element matrix
   *
   * Throws an assertion, in case any specification is missing for the type of
   * cell.
   */
  ElemMat Eval(const lf::mesh::Entity& cell) override;

  [[nodiscard]] size_type TrialComponent() const override {
    return trial_component_;
  }

  [[nodiscard]] size_type TestComponent() const override {
    return test_component_;
  }

  ~ConvectionElementMatrixProvider() override = default;

 private:
  /** @brief functor providing "first" convection coefficient */
  CONVECTION_COEFF_1 beta_1_;
  /** @brief functor providing the "second" convection coefficient */
  CONVECTION_COEFF_2 beta_2_;

  /** @brief array containing the precomputed information for the trial space:
   * fe_precomp_trial[i] contains the precomputed reference finite element of
   * the trial space component u for ref_el i.*/
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_trial_;
  /** @brief array containing the precomputed information for the test space:
   * fe_precomp_trial[i] contains the precomputed reference finite element of
   * the test space component v for ref_el i.*/
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_test_;

  /** @brief index of the involved trial component */
  size_type trial_component_;
  /** @brief index of the involved test component */
  size_type test_component_;
};

// template deduction hint.
template <class PTR, typename CONVECTION_COEFF_1, typename CONVECTION_COEFF_2>
ConvectionElementMatrixProvider(PTR fe_trial, PTR fe_test,
                                size_type trial_component,
                                size_type test_component,
                                CONVECTION_COEFF_1 beta_1,
                                CONVECTION_COEFF_2 beta_2)
    ->ConvectionElementMatrixProvider<typename PTR::element_type::SCALAR,
                                      CONVECTION_COEFF_1, CONVECTION_COEFF_2>;

// Constructor. internal construction of quadrature rules.
template <typename SCALAR, typename CONVECTION_COEFF_1,
          typename CONVECTION_COEFF_2>
ConvectionElementMatrixProvider<SCALAR, CONVECTION_COEFF_1,
                                CONVECTION_COEFF_2>::
    ConvectionElementMatrixProvider(
        std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
        std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
        size_type trial_component, size_type test_component,
        CONVECTION_COEFF_1 beta_1, CONVECTION_COEFF_2 beta_2)
    : beta_1_(std::move(beta_1)),
      beta_2_(std::move(beta_2)),
      trial_component_(trial_component),
      test_component_(test_component),
      fe_precomp_trial_(),
      fe_precomp_test_() {
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    // obtain descriptions of local shape functions.
    auto fe_trial =
        fe_space_trial->ShapeFunctionLayout(ref_el, trial_component_);
    auto fe_test = fe_space_test->ShapeFunctionLayout(ref_el, test_component_);

    // check, that  shape functions for both the trial and test space component
    // for that entity are available.
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe_trial != nullptr && fe_test != nullptr) {
      // use a quadrature rule whose degree is the sum of the degrees of the fe
      // spaces.
      size_type degree = fe_trial->Degree() + fe_test->Degree();
      lf::quad::QuadRule qr = lf::quad::make_QuadRule(ref_el, degree);
      // Precompute cell-independent quantities.
      fe_precomp_trial_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(fe_trial, qr);
      fe_precomp_test_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(fe_test, qr);
    }
  }
}

template <typename SCALAR, typename CONVECTION_COEFF_1,
          typename CONVECTION_COEFF_2>
typename ConvectionElementMatrixProvider<SCALAR, CONVECTION_COEFF_1,
                                         CONVECTION_COEFF_2>::ElemMat
ConvectionElementMatrixProvider<
    SCALAR, CONVECTION_COEFF_1,
    CONVECTION_COEFF_2>::Eval(const lf::mesh::Entity& cell) {
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>& pfe_trial =
      fe_precomp_trial_[cell.RefEl().Id()];
  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>& pfe_test =
      fe_precomp_test_[cell.RefEl().Id()];

  // check initialization:
  LF_ASSERT_MSG(pfe_trial.isInitialized(),
                "No local shape function on trial space for entity type "
                    << cell.RefEl());
  LF_ASSERT_MSG(
      pfe_test.isInitialized(),
      "No local shape function on test space for entity type " << cell.RefEl());
  // check the consistency of the quadratrue rules used.
  LF_ASSERT_MSG(pfe_trial.Qr().Points() == pfe_test.Qr().Points() &&
                    pfe_trial.Qr().Weights() == pfe_test.Qr().Weights() &&
                    pfe_trial.Qr().NumPoints() == pfe_test.Qr().NumPoints(),
                "Qr missmatch");

  // query the geometry of the cell.
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry");
  LF_ASSERT_MSG(geo_ptr->DimLocal() == 2, "Only 2D implementation available");

  // phsical dimension of the cell.
  const lf::uscalfe::dim_t world_dim = geo_ptr->DimGlobal();

  // retrive information for parametric fem computation:
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe_trial.Qr().Points()));
  LF_ASSERT_MSG(determinants.size() == pfe_trial.Qr().NumPoints(),
                "Mismatch " << determinants.size() << " <-> "
                            << pfe_trial.Qr().NumPoints());

  const Eigen::MatrixXd JinvT(
      geo_ptr->JacobianInverseGramian(pfe_trial.Qr().Points()));
  LF_ASSERT_MSG(
      JinvT.cols() == 2 * pfe_trial.Qr().NumPoints() &&
          JinvT.rows() == world_dim,
      "Mismatch " << JinvT.cols() << " <-> " << 2 * pfe_trial.Qr().NumPoints());

  // evaluate coefficients in the quadrature points
  auto beta_1_eval = beta_1_(cell, pfe_trial.Qr().Points());
  auto beta_2_eval = beta_2_(cell, pfe_trial.Qr().Points());

  // initialize element matrix
  elem_mat_t mat(pfe_test.NumRefShapeFunctions(),
                 pfe_trial.NumRefShapeFunctions());
  mat.setZero();

  // loop over quadrature points:
  for (size_type k = 0; k < pfe_trial.Qr().NumPoints(); ++k) {
    // transformed quadrature wiehgt.
    const double w = pfe_trial.Qr().Weights()[k] * determinants[k];

    // transformed gradients.
    const auto trf_trial_grad(
        JinvT.block(0, 2 * k, world_dim, 2) *
        pfe_trial.PrecompGradientsReferenceShapeFunctions()
            .block(0, 2 * k, mat.cols(), 2)
            .transpose());
    const auto trf_test_grad(JinvT.block(0, 2 * k, world_dim, 2) *
                             pfe_test.PrecompGradientsReferenceShapeFunctions()
                                 .block(0, 2 * k, mat.rows(), 2)
                                 .transpose());
    // evaluated local shape function.
    const auto trf_trial_eval(
        pfe_trial.PrecompReferenceShapeFunctions().col(k));
    const auto trf_test_eval(pfe_test.PrecompReferenceShapeFunctions().col(k));

    // update element matrix.
    const auto beta_trf_test_grad = beta_1_eval[k].transpose() * trf_test_grad;
    const auto beta_trf_trial_grad =
        beta_2_eval[k].transpose() * trf_trial_grad;
    mat += w * (beta_trf_test_grad.transpose() * trf_trial_eval.transpose() +
                trf_test_eval * beta_trf_trial_grad);
  }
  return mat;
}

/**
 * @headerfile projects/dpg/loc_comp_dpg.h
 * @brief class for local quadrature based computations of sub element matrices
 * corresponding to element matrices associated to a flux variable in the trial
 space
 *
 *
 * @tparam SCALAR type for the entries of the element matrices. Must be a field
 *                     type such as `double` or `std::complex<double>`
 * @tparam DIFF_COEFF a \ref mesh_function "MeshFunction" that defines the
 *                    coefficient \f$ \alpha \f$. It should scalar valued.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * The element matrix corresponds to the (local) bilinear form
 * @f[
    (\hat{q}_n,v) \mapsto\int\limits_{\partial K} \hat{q}_n
 \mathrm{sgn}_K(\mathbf{x}) \alpha(\mathbf{x}) v \mathrm{d}\mathbf{x}
 \;,
 * @f]
 *
 *\f$ \hat{q_n} \f$ is a flux component of the trial space and \f$ v\f$ is a
 component of the test space. *We assume that the flux is represented by
 discontinuous (interior) local shape functions on segments.
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - DIFF_COEFF must provide a scalar-valued evaluation operator
 * `operator (const Entity &,ref_coord_t)`.
 *
 */
template <typename SCALAR, typename DIFF_COEFF>
class FluxElementMatrixProvider : public SubElementMatrixProvider<SCALAR> {
  static_assert(lf::uscalfe::isMeshFunction<DIFF_COEFF>);

 public:
  /** @brief inherited types for element matrices */
  using elem_mat_t = typename SubElementMatrixProvider<SCALAR>::elem_mat_t;
  using ElemMat = typename SubElementMatrixProvider<SCALAR>::ElemMat;

  FluxElementMatrixProvider(const FluxElementMatrixProvider&) = delete;
  FluxElementMatrixProvider(FluxElementMatrixProvider&&) noexcept = default;
  FluxElementMatrixProvider& operator=(const FluxElementMatrixProvider&) =
      delete;
  FluxElementMatrixProvider& operator=(FluxElementMatrixProvider&&) = delete;

  /**
   * @brief Constructor: performs cell-independent precomputations.
   * @param fe_space_trial collection of specifications for the trial fe space
   * @param fe_space_test collection of specification for the test fe space
   * @param trial_component Index of the trial space component \f$\hat{q}_n\f$.
   *  This component should represent a flux and
   * be described via local discontinuous (interior) shape functions on the
   * reference line segment.
   * @param test_component Index of the  test space component \f$v\f$
   * @param alpha  mesh function for the scalar valued coefficient involved in
   * the bilinear form.
   *
   * This constructur uses local quadrature rules. On each edge of the boundary
   * a quadrature rule with degree of exactness chosen as the sum of the
   * polynomial degree of the two components \f$ \hat{q}_n \f$ and \f$ v \f$ is
   * used.
   */
  FluxElementMatrixProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      size_type trial_component, size_type test_component, DIFF_COEFF alpha);

  /**
   * @brief main routine for the computation of element matrices
   * @param cell reference to a (triangular or quadrilateral) cell for which the
   *        element matirx should be computed.
   * @return small dense matrix containing the element matrix
   *
   * Actual computation is based on numerical quadrature and mapping techniques.
   *
   * Throws an assertion, in case any specification is missing for the type of
   * cell.
   */
  ElemMat Eval(const lf::mesh::Entity& cell) override;

  [[nodiscard]] size_type TrialComponent() const override {
    return trial_component_;
  }

  [[nodiscard]] size_type TestComponent() const override {
    return test_component_;
  }

  ~FluxElementMatrixProvider() override = default;

 private:
  /** @brief functor providing the coefficient invoved in the bilinear form */
  DIFF_COEFF alpha_;
  /** @brief array containing the precomputed information for the trial space
   * flux component: fe_precomp_trial[i] contains the precomputed reference
   * finite element of flux which is used in the computations
   * of element matrices on cells of ref_el i. Note that the actual precomputed
   * reference finite element is associated with ref_el
   * lf::base::RefEl::kSegment().
   */
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_trial_;
  /** @brief array containing the precomputed information for the trial space:
   * fe_precomp_trial[i][j] contains the precomputed reference finite element of
   * \f$ v \f$ for
   * ref_el i, evaluated on the quadrature points transformed to edge j of the
   * boundary of the reference element*/
  std::array<
      std::vector<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>>,
      5>
      fe_precomp_test_;

  /** @brief index of the trial component \f$ \hat{q}_n \f$ */
  size_type trial_component_;
  /** @brief index of test component v */
  size_type test_component_;
  /** @brief helper class, to evaluate the  sgn_K function*/
  PrescribedSignProvider sign_provider_;
};

// template deduction hint.
template <class PTR, class DIFF_COEFF>
FluxElementMatrixProvider(PTR fe_space_trial, PTR fe_space_test,
                          size_type trial_component, size_type test_component,
                          DIFF_COEFF alpha)
    ->FluxElementMatrixProvider<typename PTR::element_type::SCALAR, DIFF_COEFF>;

// constructor internal construction of quadrature rules.
template <typename SCALAR, typename DIFF_COEFF>
FluxElementMatrixProvider<SCALAR, DIFF_COEFF>::FluxElementMatrixProvider(
    std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
    std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
    size_type trial_component, size_type test_component, DIFF_COEFF alpha)
    : alpha_(std::move(alpha)),
      trial_component_(trial_component),
      test_component_(test_component),
      sign_provider_(fe_space_trial->Mesh()),
      fe_precomp_test_(),
      fe_precomp_trial_() {
  // obtain description of the shape functions representing the flux.
  auto fe_trial = fe_space_trial->ShapeFunctionLayout(
      lf::base::RefEl::kSegment(), trial_component_);
  LF_ASSERT_MSG(fe_trial != nullptr, "Missing description of flux space");
  LF_ASSERT_MSG(fe_trial->NumRefShapeFunctions(1) == 0,
                "shape functions associated with endpoints not supported.");
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    // obtain description of the local shape functions in the test space.
    auto fe_test = fe_space_test->ShapeFunctionLayout(ref_el, test_component_);

    // check, that  shape functions of test space component for
    // that entity are available.
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe_test != nullptr) {
      // use a quadrature rule whose degree is the sum of the degrees of the fe
      // spaces.
      int degree = fe_trial->Degree() + fe_test->Degree();
      lf::quad::QuadRule segment_qr =
          lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), degree);

      // precompute cell-independent quantities for the flux
      fe_precomp_trial_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(fe_trial,
                                                               segment_qr);

      // transform quadrule to all parts of the boundary.
      std::vector<lf::quad::QuadRule> boundaryQr =
          BoundaryQuadRule(ref_el, segment_qr);
      int numSegments = ref_el.NumSubEntities(1);
      LF_ASSERT_MSG(boundaryQr.size() == numSegments,
                    "boundaryQr.size() = " << boundaryQr.size() << " <-> "
                                           << "numSegments= " << numSegments);
      fe_precomp_test_[ref_el.Id()].resize(numSegments);

      for (int segment = 0; segment < numSegments; segment++) {
        // On each edge of the boundary, compute cell-independent quantities for
        // the test soace,
        fe_precomp_test_[ref_el.Id()][segment] =
            lf::uscalfe::PrecomputedScalarReferenceFiniteElement(
                fe_test, boundaryQr[segment]);
      }
    }
  }
}

template <typename SCALAR, typename DIFF_COEFF>
typename FluxElementMatrixProvider<SCALAR, DIFF_COEFF>::ElemMat
FluxElementMatrixProvider<SCALAR, DIFF_COEFF>::Eval(
    const lf::mesh::Entity& cell) {
  // ibtain precomputed information
  auto& pfe_trial = fe_precomp_trial_[cell.RefEl().Id()];
  auto& pfes_test = fe_precomp_test_[cell.RefEl().Id()];
  // the number of edges of the cell boundary
  int numSegments = cell.RefEl().NumSubEntities(1);

  // check initialization:
  LF_ASSERT_MSG(pfes_test.size() == numSegments,
                "Missing local shape functions for test space on ref_el "
                    << cell.RefEl());
  LF_ASSERT_MSG(pfe_trial.isInitialized(),
                "missing local shape functions for trial space on ref_el"
                    << cell.RefEl());
  for (int i = 0; i < numSegments; ++i) {
    LF_ASSERT_MSG(pfes_test[i].isInitialized(),
                  "missing local shape functions for test space on ref_el "
                      << cell.RefEl());
    LF_ASSERT_MSG(pfes_test[i].Qr().NumPoints() == pfe_trial.Qr().NumPoints(),
                  "qr missmatch between trial and test space.");
  }

  // query the geometry of the cell.
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry");
  LF_ASSERT_MSG(geo_ptr->DimLocal() == 2, "Only 2D implementation available");

  // initialize element matrix
  elem_mat_t mat(pfes_test[0].NumRefShapeFunctions(),
                 pfe_trial.NumRefShapeFunctions() * numSegments);
  mat.setZero();

  // integrate over all parts of the boundary,
  for (int segment = 0; segment < numSegments; segment++) {
    // query the geometry of the segment.
    const auto segment_geo_ptr = geo_ptr->SubGeometry(1, segment);
    LF_ASSERT_MSG(segment_geo_ptr != nullptr, "Invalid geometry");

    // retrive determinants for paramteric fem
    const Eigen::VectorXd determinants =
        segment_geo_ptr->IntegrationElement(pfe_trial.Qr().Points());
    LF_ASSERT_MSG(determinants.size() == pfe_trial.Qr().NumPoints(),
                  "Mismatch " << determinants.size() << " <-> "
                              << pfe_trial.Qr().NumPoints());

    // evaluate coefficient- and sign function.
    auto alphaeval = alpha_(cell, pfes_test[segment].Qr().Points());
    int sgn =
        sign_provider_.PrescribedSign(cell, *cell.SubEntities(1)[segment]);

    // iterate over quadrature points.
    for (size_type k = 0; k < pfe_trial.Qr().NumPoints(); ++k) {
      // transformed quadrature weight
      const double w = pfe_trial.Qr().Weights()[k] * determinants[k] * sgn;
      const auto trial_eval = pfe_trial.PrecompReferenceShapeFunctions().col(k);
      const auto test_eval =
          pfes_test[segment].PrecompReferenceShapeFunctions().col(k);

      // use the assumptions about the local shape function layout of the flux
      // component and the rules for ordering local shape functions, to update
      // correct part of the element matrix.
      mat.block(0, segment * pfe_trial.NumRefShapeFunctions(),
                pfes_test[segment].NumRefShapeFunctions(),
                pfe_trial.NumRefShapeFunctions()) +=
          alphaeval[k] * w * test_eval * trial_eval.transpose();
    }
  }
  return mat;
}

/**
 *
 * @headerfile projects/dpg/loc_comp_dpg.h
 * @brief Class for local quadrature based computations of sub element matrices
 corresponding
 *  to element matrices associated to a trace variable in the trial space
 *
 * @tparam SCALAR type for the entries of the element matrices. Must be a field
 *                     type such as `double` or `std::complex<double>`
 * @tparam COEFF a \ref mesh_function "MeshFunction" that defines a coefficient
 * \f$ \mathbf{\beta} \f$. It should be vector valued.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * The element matrix corresponds to the (local) bilinear form
 * @f[
    (\hat u,v) \mapsto\int\limits_{\partial K} \hat u \mathbf{n}_K \cdot
 \boldsymbol{\beta}(\mathbf{x}) v\mathrm{d}\mathbf{x}
 \;,
 * @f]
 * \f$ \hat{u} \f$  is a trace  component of the trial space and \f$ v \f$ a
 component of the test space.

 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - COEFF must provide an evaluation operator
 * `operator (const Entity &,ref_coord_t)` returning a vector, that is
 compatible with Eigen's matrices. Usually it will
 * be an Eigen::Vector of variable or fixed size.
 *
 */
template <typename SCALAR, typename COEFF>
class TraceElementMatrixProvider : public SubElementMatrixProvider<SCALAR> {
  static_assert(lf::uscalfe::isMeshFunction<COEFF>);

 public:
  /** @brief inherited types for element matrices */
  using elem_mat_t = typename SubElementMatrixProvider<SCALAR>::elem_mat_t;
  using ElemMat = typename SubElementMatrixProvider<SCALAR>::ElemMat;

  TraceElementMatrixProvider(const TraceElementMatrixProvider&) = delete;
  TraceElementMatrixProvider(TraceElementMatrixProvider&&) noexcept = default;
  TraceElementMatrixProvider& operator=(const TraceElementMatrixProvider&) =
      delete;
  TraceElementMatrixProvider& operator=(TraceElementMatrixProvider&&) = delete;

  /**
   * @brief Constructor: performs cell-independent precomputations.
   * @param fe_space_trial collection of specifications for the trial fe space
   * @param fe_space_test collection of specifications for the test fe space
   * @param trial_component Index of the trial space component \f$\hat u\f$
   * @param test_component Index of the  test space component \f$v\f$
   * @param beta mesh function for the vector valued coefficient involved in the
   * bilinear form.
   *
   * This constructur uses local quadrature rules. On each edge of the boundary
   * a quadrature rule with degree of exactness chosen as the sum of the
   * polynomial degree of the two components \f$ \hat{u} \f$ and \f$ v \f$ is
   * used.
   */
  TraceElementMatrixProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      size_type trial_component, size_type test_component, COEFF beta);
  /**
   * @brief main routine for the computation of element matrices
   * @param cell reference to a (triangular or quadrilateral) cell for which the
   *        element matirx should be computed.
   * @return small dense matrix containing the element matrix
   *
   * Actual computation is based on numerical quadrature and mapping techniques.
   *
   * Throws an assertion, in case any specification is missing for the type of
   * cell.
   */
  ElemMat Eval(const lf::mesh::Entity& cell) override;

  [[nodiscard]] size_type TrialComponent() const override {
    return trial_component_;
  }

  [[nodiscard]] size_type TestComponent() const override {
    return test_component_;
  }

  ~TraceElementMatrixProvider() override = default;

 private:
  /** \brief functor providing the coefficient invoved in the bilinear form */
  COEFF beta_;

  /** @brief array containing the precomputed information for the trial space:
   * fe_precomp_trial[i][j] contains the precomputed reference finite element of
   * \f$ \hat{u} \f$  for
   * ref_el i  evaluated on the quadrature points transformed to edge j of the
   * boundary.*/
  std::array<
      std::vector<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>>,
      5>
      fe_precomp_trial_;
  /** @brief array containing the precomputed information for the test space:
   * fe_precomp_trial[i][j] contains the precomputed reference finite element of
   * \f$ v \f$  for
   * ref_el i evaluated on the quadrature points transformed to edge j of the
   * boundary.*/
  std::array<
      std::vector<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>>,
      5>
      fe_precomp_test_;

  /** @brief the original segment quad rule used to construct the quadrules on
   * ref_el boundaries. segment_qr[i] contains the segment qr for ref_el i. */
  std::array<lf::quad::QuadRule, 5> segment_qr_;

  /** @brief index of trial component \f$ \hat{u} \f$ */
  size_type trial_component_;
  /** @brief index of test component v */
  size_type test_component_;
};

// template deduction hint.
template <class PTR, class COEFF>
TraceElementMatrixProvider(PTR fe_space_trial, PTR fe_space_test,
                           size_type trial_component, size_type test_component,
                           COEFF alpha)
    ->TraceElementMatrixProvider<typename PTR::element_type::SCALAR, COEFF>;

// Constructor. internal construction of quadrature rules.
template <typename SCALAR, typename COEFF>
TraceElementMatrixProvider<SCALAR, COEFF>::TraceElementMatrixProvider(
    std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_trial,
    std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
    size_type trial_component, size_type test_component, COEFF beta)
    : beta_(std::move(beta)),
      trial_component_(trial_component),
      test_component_(test_component),
      fe_precomp_test_(),
      fe_precomp_trial_() {
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    // obtain descriptions of local shape functions.
    auto fe_trial =
        fe_space_trial->ShapeFunctionLayout(ref_el, trial_component_);
    auto fe_test = fe_space_test->ShapeFunctionLayout(ref_el, test_component_);

    // check, that  shape functions for both the trial and test space component
    // for that entity are available.
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe_trial != nullptr && fe_test != nullptr) {
      // use a quadrature rule whose degree is the sum of the degrees of the fe
      // spaces.
      size_type degree = fe_trial->Degree() + fe_test->Degree();
      // initialize a quadrule of a the given degree on the reference line
      // segment.
      lf::quad::QuadRule segment_qr =
          lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), degree);
      segment_qr_[ref_el.Id()] = segment_qr;

      // transform it to all parts of the boundary.
      std::vector<lf::quad::QuadRule> boundary_qr =
          BoundaryQuadRule(ref_el, segment_qr);
      int num_segments = ref_el.NumSubEntities(1);
      LF_ASSERT_MSG(
          boundary_qr.size() == num_segments,
          "boundary_qr.size() = " << boundary_qr.size() << "<->"
                                  << "num_segments = " << num_segments);

      fe_precomp_trial_[ref_el.Id()].resize(num_segments);
      fe_precomp_test_[ref_el.Id()].resize(num_segments);

      for (int segment = 0; segment < num_segments; segment++) {
        // on each edge of the boundary precompute cell-independent quantities.
        fe_precomp_trial_[ref_el.Id()][segment] =
            lf::uscalfe::PrecomputedScalarReferenceFiniteElement(
                fe_trial, boundary_qr[segment]);
        fe_precomp_test_[ref_el.Id()][segment] =
            lf::uscalfe::PrecomputedScalarReferenceFiniteElement(
                fe_test, boundary_qr[segment]);
      }
    }
  }
}

template <typename SCALAR, typename COEFF>
typename TraceElementMatrixProvider<SCALAR, COEFF>::ElemMat
TraceElementMatrixProvider<SCALAR, COEFF>::Eval(const lf::mesh::Entity& cell) {
  // obtain precomputed about values of local shape functions
  // and their gradients at quadrature points on all parts of the boundary.
  auto& pfes_trial = fe_precomp_trial_[cell.RefEl().Id()];
  auto& pfes_test = fe_precomp_test_[cell.RefEl().Id()];

  // retrive the original quadrule on the segment, used only for the weights.
  auto& segment_qr = segment_qr_[cell.RefEl().Id()];

  // the number of edges of the cell boundary.
  size_type num_segments = cell.RefEl().NumSubEntities(1);

  // check initialization:
  LF_ASSERT_MSG(pfes_trial.size() == num_segments,
                "Missing local shape functions on trial space for ref_el"
                    << cell.RefEl());
  LF_ASSERT_MSG(
      pfes_test.size() == num_segments,
      "Missing local shape functions on test space for ref_el" << cell.RefEl());

  // check consistency of quadrqture rules between fe spaces.
  for (size_type segment = 0; segment < num_segments; segment++) {
    LF_ASSERT_MSG(pfes_trial[segment].isInitialized(),
                  "missing local shape functions on trial space for ref_el"
                      << cell.RefEl());
    LF_ASSERT_MSG(pfes_test[segment].isInitialized(),
                  "missing local shape functions on trial space for ref_el"
                      << cell.RefEl());
    LF_ASSERT_MSG(
        pfes_trial[segment].Qr().Points() == pfes_test[segment].Qr().Points() &&
            pfes_trial[segment].Qr().Weights() ==
                pfes_test[segment].Qr().Weights() &&
            pfes_trial[segment].Qr().NumPoints() ==
                pfes_test[segment].Qr().NumPoints(),
        "Qr missmatch");
  }

  // check consistency with semgment quad rule:
  LF_ASSERT_MSG(segment_qr.NumPoints() == pfes_trial[0].Qr().NumPoints(),
                "QR missmatch");

  // query the geometry of the cell.
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry");
  LF_ASSERT_MSG(geo_ptr->DimLocal() == 2, "Only 2D implementation available");

  // initialize element matrix:
  elem_mat_t mat(pfes_test[0].NumRefShapeFunctions(),
                 pfes_trial[0].NumRefShapeFunctions());
  mat.setZero();

  // retrive normals on all edges of the cell
  Eigen::MatrixXd normals = OuterNormals(*geo_ptr);
  LF_ASSERT_MSG(normals.cols() == num_segments,
                "normals.cols()= " << normals.cols() << "<->"
                                   << "num_segments= " << num_segments);

  // sum up by integrating over all parts of the boundary.
  for (size_type segment = 0; segment < num_segments; segment++) {
    // query the geometry of the segment
    const auto segment_ptr = geo_ptr->SubGeometry(1, segment);
    LF_ASSERT_MSG(segment_ptr != nullptr, "Invalid geometry");

    // retrive precomputed information on current edge.
    auto& pfe_trial = pfes_trial[segment];
    auto& pfe_test = pfes_test[segment];

    // retrive determinants for paramteric fem.
    const Eigen::VectorXd determinants(
        segment_ptr->IntegrationElement(segment_qr.Points()));
    LF_ASSERT_MSG(determinants.size() == pfe_trial.Qr().NumPoints(),
                  "Mismatch " << determinants.size() << " <-> "
                              << segment_qr.NumPoints());

    // evaluate coefficient function
    auto betaeval = beta_(cell, pfe_trial.Qr().Points());

    // iterate over quadrule points on current edge of the boundary.
    for (size_type k = 0; k < pfe_trial.Qr().NumPoints(); ++k) {
      // transformed quadrature weight
      const double w = segment_qr.Weights()[k] * determinants[k];

      const auto trial_eval = pfe_trial.PrecompReferenceShapeFunctions().col(k);
      const auto test_eval = pfe_test.PrecompReferenceShapeFunctions().col(k);
      const auto beta_dot_n = normals.col(segment).transpose() * betaeval[k];

      // update element matrix.
      mat += w * test_eval * trial_eval.transpose() * beta_dot_n;
    }
  }
  return mat;
}

/**
 * @headerfile projects/dpg/loc_comp_dpg.h
 * @brief Class for local quadrature based computations of sub vectors
 corresponding to
 * load vectors.
 *
 * @tparam SCALAR type for the entries of the element vectors. Must be a
 * fiel type such as 'double'.
 * @tparam FUNCTOR a \ref mesh_function "MeshFunction" that defines the source
 * function \f$ f \f$. It should be scalar-valued.
 *
 * @note This class complies with the type requirements for the template
 argument
 * 'ELEM_VEC_COMP' of the function
 * lf::assemble::AssembleVectorLocally89
 *
 *
 * The load vector corresponds to the local linear form
 * @f[
      v \mapsto \int_K f(\mathbf{x})\,v\,\mathrm{d}\mathbf{x}\;,
 * @f]
 * with source function \f$ f \f$. \f$v \f$ is a component of the test space.
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like double
 * - FUNCTOR must provide an evaluation operator
 * `operator (const Entity &,ref_coord_t)` that returns a scalar
 */
template <typename SCALAR, typename FUNCTOR>
class LoadElementVectorProvider : public SubElementVectorProvider<SCALAR> {
  static_assert(lf::uscalfe::isMeshFunction<FUNCTOR>);

 public:
  /** @brief inherited types for element vectors */
  using elem_vec_t = typename SubElementVectorProvider<SCALAR>::elem_vec_t;
  using ElemVec = typename SubElementVectorProvider<SCALAR>::ElemVec;

  LoadElementVectorProvider(const LoadElementVectorProvider&) = delete;
  LoadElementVectorProvider(LoadElementVectorProvider&&) noexcept = default;
  LoadElementVectorProvider& operator=(const LoadElementVectorProvider&) =
      delete;
  LoadElementVectorProvider& operator=(LoadElementVectorProvider&&) = delete;

  /**
   * @brief Constructor: performs cell-independent precomputations.
   * @param fe_space_test collection of specifications for the fe space
   * @param test_component Index of the  test space component \f$v\f$
   * @param f mesh function for the scalar valued source function
   *
   * This constructor uses local quadrature rules with twice the degree of
   * exactnes of the polynomial degree of \f$ v\f$
   */
  LoadElementVectorProvider(
      std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
      size_type test_component, FUNCTOR f);

  /**
   * @brief main routine for the computation of element vectors
   * @param cell reference to a (triangular or quadrilateral) cell for which the
   * element vector should be computed
   * @return small column vector containing the element vector.
   *
   */
  ElemVec Eval(const lf::mesh::Entity& cell) override;

  [[nodiscard]] size_type TestComponent() const override {
    return test_component_;
  }

  ~LoadElementVectorProvider() override = default;

 private:
  /** @brief functor providing the source function */
  FUNCTOR f_;
  /** @brief array containing precomputed information for the test space:
   * fe_precomp_test_[i] contains the precomputed reference finite element for
   * \f$ v\f$ on ref_el i. */
  std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<SCALAR>, 5>
      fe_precomp_test_;

  /** @brief index of test component v */
  size_type test_component_;
};

// template deduction hint
template <class PTR, class FUNCTOR>
LoadElementVectorProvider(PTR fe_space_test, size_type test_component,
                          FUNCTOR f)
    ->LoadElementVectorProvider<typename PTR::element_type::SCALAR, FUNCTOR>;

// constructor: internal construction of quadrature rules.
template <typename SCALAR, typename FUNCTOR>
LoadElementVectorProvider<SCALAR, FUNCTOR>::LoadElementVectorProvider(
    std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space_test,
    size_type test_component, FUNCTOR f)
    : f_(std::move(f)), test_component_(test_component), fe_precomp_test_() {
  for (auto ref_el : {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}) {
    // obtain description of local shape function
    auto fe_test = fe_space_test->ShapeFunctionLayout(ref_el, test_component_);
    // check, that  shape functions for this component are available
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe_test != nullptr) {
      // Precompute cell-independent quantities based on quadrature rules
      // with twice the degree of exactness compared to the degree of the
      // finite element space.
      fe_precomp_test_[ref_el.Id()] =
          lf::uscalfe::PrecomputedScalarReferenceFiniteElement(
              fe_test, lf::quad::make_QuadRule(ref_el, 2 * fe_test->Degree()));
    }
  }
}

template <typename SCALAR, typename FUNCTOR>
typename LoadElementVectorProvider<SCALAR, FUNCTOR>::ElemVec
LoadElementVectorProvider<SCALAR, FUNCTOR>::Eval(const lf::mesh::Entity& cell) {
  // Obtain precomputed information about values of local shape functions and
  // their gradients at quadrature points.
  auto& pfe_test = fe_precomp_test_[cell.RefEl().Id()];

  // check initilization:
  LF_ASSERT_MSG(pfe_test.isInitialized(),
                "missing local shape function on ref_el " << cell.RefEl());

  // querry geometry
  const lf::geometry::Geometry* geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry");
  LF_ASSERT_MSG(geo_ptr->DimLocal() == 2, "Only 2D implementation available");

  // retrive information for parametric fem computation
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe_test.Qr().Points()));
  LF_ASSERT_MSG(determinants.size() == pfe_test.Qr().NumPoints(),
                "Mismatch " << determinants.size() << " <-> "
                            << pfe_test.Qr().NumPoints());

  // evaluate f at the quadrature points
  auto fval = f_(cell, pfe_test.Qr().Points());

  // initialize the element vector
  elem_vec_t vec(pfe_test.NumRefShapeFunctions());
  vec.setZero();

  // loop over quadrature points:
  for (size_type k = 0; k < pfe_test.Qr().NumPoints(); ++k) {
    // transformed quadrature weight
    const double w = pfe_test.Qr().Weights()[k] * determinants[k];
    vec += fval[k] * w * pfe_test.PrecompReferenceShapeFunctions().col(k);
  }

  return vec;
}

}  // namespace projects::dpg
#endif  // PROJECTS_DPG_LOC_COMP_DPG

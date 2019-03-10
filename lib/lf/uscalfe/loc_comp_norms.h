#ifndef LF_LOCCOMPNORMS_H
#define LF_LOCCOMPNORMS_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Classes for the local computation of norms of (differences of)
 * (FE) functions
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <cmath>
#include <iostream>
#include <type_traits>

#include <lf/quad/quad.h>
#include "loc_comp_ellbvp.h"

namespace lf::uscalfe {
/**
 * @brief Class for the local computation of the L2 norm of the difference of
 *        a Lagrangian finite element function and a generic function
 *
 * @tparam FUNCTOR describes the function from which the finite element function
 * is subtracted, should model the concept \ref mesh_function .
 *
 * This class complies with the type requirements for the template argument
 * LOCAL_NORM_COMP of the function lf::uscalfe::L2NormDifference.
 * Its main method is the evaluation operator `operator ()`.
 *
 * ### Use case
 *
 * ~~~
 *   MeshFunctionL2NormDifference loc_comp(
 *     fe_space,[](auto x) { return (x[0] * (1.0 - x[1])); }, quad_order);
 * ~~~
 */
template <typename FUNCTOR>
class MeshFunctionL2NormDifference {
  static_assert(isMeshFunction<FUNCTOR>);

 public:
  /** @brief standard constructors */
  /** @{ */
  MeshFunctionL2NormDifference(const MeshFunctionL2NormDifference &) = delete;
  MeshFunctionL2NormDifference(MeshFunctionL2NormDifference &&) noexcept =
      default;
  MeshFunctionL2NormDifference &operator=(
      const MeshFunctionL2NormDifference &) = delete;
  MeshFunctionL2NormDifference &operator=(MeshFunctionL2NormDifference &&) =
      delete;
  /** @} */

  /** @brief Constructor
   *
   * @param fe_space Describes the Shapefunctions over which to integrate.
   * @param u object providing scalar-valued function
   * @param loc_quad_degree desired degree of exctness of local quadrature,
   * default value = 0. If = 0, the quadrature degree is determined from the
   * polynomial degree of the reference shape functions.
   *
   * @note The finite element space need not provide reference shape functions
   *       for all cell types if these cell types do not occur.
   */
  MeshFunctionL2NormDifference(
      const std::shared_ptr<const UniformScalarFESpace<double>> &fe_space,
      FUNCTOR u, quad::quadDegree_t loc_quad_order)
      : u_(std::move(u)) {
    LF_ASSERT_MSG(fe_space != nullptr, "Invalid FE space passed");
    for (auto ref_el : {base::RefEl::kTria(), base::RefEl::kQuad()}) {
      auto fe = fe_space->ShapeFunctionLayout(ref_el);
      if (fe != nullptr) {
        fe_precomp_[ref_el.Id()] =
            PrecomputedScalarReferenceFiniteElement<double>(
                std ::move(fe), quad ::make_QuadRule(ref_el, loc_quad_order));
      }
    }
  }

  /**
   * @brief main function for the computation of **squared** local L2 norms
   *
   * @tparam DOFVECTOR vector type for passing local coefficient vector
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @param dofs vector of degrees of freedom for the current entity
   * @return the local "error norm" squared
   *
   * Actual computation of the element matrix based on numerical quadrature and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying Lagrangian finite element spaces: for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials o degree 2p.
   *
   * ### Template parameter type requirements
   *
   * - DOFVECTOR must be compatible with std::vector<SCALAR>
   */
  template <typename DOFVECTOR>
  double operator()(const lf::mesh::Entity &cell, const DOFVECTOR &dofs);

  /** Virtual destructor */
  virtual ~MeshFunctionL2NormDifference() = default;

 private:
  /** @brief a handle to the generic function with which to compare the finite
   * element function */
  FUNCTOR u_;

  /** @brief cell-type dependent but cell-independent information */
  std::array<PrecomputedScalarReferenceFiniteElement<double>, 5> fe_precomp_;

 public:
  /** @brief output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_cell = 8;
  static const unsigned int kout_comp = 16;
};

// Initialization of control variable
template <typename FUNCTOR>
unsigned int MeshFunctionL2NormDifference<FUNCTOR>::ctrl_ = 0;

template <typename FUNCTOR>
template <typename DOFVECTOR>
double MeshFunctionL2NormDifference<FUNCTOR>::operator()(
    const lf::mesh::Entity &cell, const DOFVECTOR &dofs) {
  using dof_t = typename DOFVECTOR::value_type;
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  auto &fe = fe_precomp_[ref_el.Id()];
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  LF_ASSERT_MSG(dofs.size() >= fe.NumRefShapeFunctions(),
                "Dof vector too short");
  // Obtain cell--type depend invariant data
  const Eigen::MatrixXd &qr_Points{fe.Qr().Points()};
  const Eigen::VectorXd &qr_Weights{fe.Qr().Weights()};
  const size_type qr_NumPts = fe.Qr().NumPoints();
  const size_type num_LSF = fe.NumRefShapeFunctions();
  const auto &rsf_QuadPts{fe.PrecompReferenceShapeFunctions()};

  SWITCHEDSTATEMENT(ctrl_, kout_cell,
                    std::cout << ref_el << ", (Nlsf = " << num_LSF
                              << ", Nqr = " << qr_NumPts << ") shape = \n"
                              << geo_ptr->Global(ref_el.NodeCoords())
                              << std::endl);

  // World coordinates of quadrature points
  const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(qr_Points));
  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr_Points));
  SWITCHEDSTATEMENT(
      ctrl_, kout_comp, std::cout << "\t dofs = "; for (auto x
                                                        : dofs) {
        std::cout << x << " ";
      } std::cout << "\t dets "
                  << determinants.transpose() << std::endl);

  double sum = 0.0;  // summation variable

  auto uvals = u_(cell, qr_Points);

  // Loop over quadrature points
  for (int k = 0; k < qr_NumPts; ++k) {
    // Value of finite element function at quadrature point
    // Loop over local shape functions/dofs to compute the value of
    // the finite element function in the quadrature point
    dof_t uh_val = 0.0;
    for (int j = 0; j < num_LSF; ++j) {
      uh_val += rsf_QuadPts(j, k) * dofs[j];
    }
    // form the difference
    uh_val -= uvals[k];
    // sum the quared modulus weighted with quadrature weight and metric factor
    sum += qr_Weights[k] * determinants[k] * std::fabs(uh_val * uh_val);
    SWITCHEDSTATEMENT(ctrl_, kout_comp,
                      std::cout << "\t @ " << (mapped_qpts.col(k)).transpose()
                                << ": uh = " << uh_val << ", u = " << uvals[k]
                                << std::endl);
  }
  SWITCHEDSTATEMENT(ctrl_, kout_comp,
                    std::cout << "\t loc norm^2 = " << sum << std::endl);
  return sum;
}

// **********************************************************************

/**
 * @ingroup mesh_function
 * @brief Class for the local computation of the L2 norm of the difference of
 *        the _gradient_ of a Lagrangian finite element function and a
 *        generic vector field
 *
 * @tparam VEC_FUNC type for a generic vector field
 *
 * ### Template parameter type requirements
 *
 * - VEC_FUNC must behave like `std::function<VECTOR(Eigen::DenseBase)>`
 * where VECTOR is a vector type like a vector of Eigen or std::vector.
 *
 * This class complies with the type requirements for the template argument
 * LOCAL_NORM_COMP of the function lf::uscalfe::L2NormDifference. Its main
 method
 * is the evaluation operator returning the **squared** local norm.
 *
 *
 * ### Use case
 *
 * ~~~
 * MeshFunctionL2GradientDifference lc_H1(fe_space,
 *    [](auto x) { return Eigen::Vector2d(-x[1],x[0]); },2)
 * ~~~
 */
template <typename VEC_FUNC>
class MeshFunctionL2GradientDifference {
  static_assert(isMeshFunction<VEC_FUNC>);

 public:
  /** defgroup stdc
   * @brief standard constructors
   * @{ */
  MeshFunctionL2GradientDifference(const MeshFunctionL2GradientDifference &) =
      delete;
  MeshFunctionL2GradientDifference(
      MeshFunctionL2GradientDifference &&) noexcept = default;
  MeshFunctionL2GradientDifference &operator=(
      const MeshFunctionL2GradientDifference &) = delete;
  MeshFunctionL2GradientDifference &operator=(
      MeshFunctionL2GradientDifference &&) = delete;
  /** @} */

  /** @brief Constructor
   *
   * @param fe_space Describes the shape functions that should be used.
   * @param vecfield object providing a generic vector field
   * @param loc_quad_order desired degree of exactness of the local quadrature
   * rule
   *
   * If loc_quad_order == 0, the required degree of exactness is deduced from
   * the polynomial degree of the finite element space.
   */
  MeshFunctionL2GradientDifference(
      const std::shared_ptr<const UniformScalarFESpace<double>> &fe_space,
      VEC_FUNC vecfield, lf::quad::quadDegree_t loc_quad_order)
      : vecfield_(std::move(vecfield)) {
    LF_ASSERT_MSG(fe_space != nullptr, "Invalid FE space passed!");
    for (auto ref_el : {base::RefEl::kTria(), base::RefEl::kQuad()}) {
      auto sfl{fe_space->ShapeFunctionLayout(ref_el)};
      if (sfl != nullptr) {
        fe_precomp_[ref_el.Id()] = PrecomputedScalarReferenceFiniteElement(
            sfl, quad::make_QuadRule(ref_el, loc_quad_order));
      }
    }
  }
  /**
   * @brief main function for the computation of **squared** local L2 norms
   *      of gradients
   *
   * @tparam DOFVECTOR vector type for passing local coefficient vector
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @param dofs vector of degrees of freedom for the current entity
   * @return the local "error norm" squared
   *
   * Actual computation of the element matrix based on numerical quadrature
   * and mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying Lagrangian finite element spaces:
   * for polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials of degree 2p ("moderate overkill quadrature")
   *
   * ### Template parameter type requirements
   *
   * - DOFVECTOR must be compatible with std::vector of a SCALAR type
   */
  template <typename DOFVECTOR>
  double operator()(const lf::mesh::Entity &cell, const DOFVECTOR &dofs);

  virtual ~MeshFunctionL2GradientDifference() = default;

 private:
  /** @brief a handle to the generic function with which to compare the
  finite
   * element function */
  VEC_FUNC vecfield_;

  std::array<PrecomputedScalarReferenceFiniteElement<double>, 5> fe_precomp_;

 public:
  /** @brief output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_cell = 8;
  static const unsigned int kout_comp = 16;
};

// Initialization of control variable
template <typename VEC_FUNC>
unsigned int MeshFunctionL2GradientDifference<VEC_FUNC>::ctrl_ = 0;

template <typename VEC_FUNC>
template <typename DOFVECTOR>
double MeshFunctionL2GradientDifference<VEC_FUNC>::operator()(
    const lf::mesh::Entity &cell, const DOFVECTOR &dofs) {
  using dof_t = typename DOFVECTOR::value_type;
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  auto &pfe = fe_precomp_[ref_el.Id()];
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  LF_ASSERT_MSG(dofs.size() >= pfe.NumRefShapeFunctions(),
                "Dof vector too short");

  // Physical dimension of the cell
  const dim_t world_dim = geo_ptr->DimGlobal();

  // Obtain cell--type depend but cell-independent data
  auto &qr_Points{pfe.Qr().Points()};
  auto &qr_Weights{pfe.Qr().Weights()};
  size_type qr_NumPts = pfe.Qr().NumPoints();
  size_type num_LSF = pfe.NumRefShapeFunctions();
  auto gradrsf_QuadPts = pfe.PrecompGradientsReferenceShapeFunctions();

  SWITCHEDSTATEMENT(ctrl_, kout_cell,
                    std::cout << ref_el << ", (Nlsf = " << num_LSF
                              << ", Nqr = " << qr_NumPts << ") shape = \n"
                              << geo_ptr->Global(ref_el.NodeCoords())
                              << std::endl);

  // World coordinates of quadrature points
  const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(qr_Points));
  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr_Points));
  // Fetch the transformation matrices for gradients
  const Eigen::MatrixXd JinvT(geo_ptr->JacobianInverseGramian(qr_Points));

  SWITCHEDSTATEMENT(
      ctrl_, kout_comp, std::cout << "\t dofs = "; for (auto x
                                                        : dofs) {
        std::cout << x << " ";
      } std::cout << "\t dets "
                  << determinants.transpose() << std::endl);

  double sum = 0.0;  // summation variable

  auto vfval = vecfield_(cell, qr_Points);

  // Loop over quadrature points
  for (int k = 0; k < qr_NumPts; ++k) {
    // Value of the gradient of the finite element function at quadrature point
    // Transformed gradients
    const auto trf_grad(
        JinvT.block(0, 2 * k, world_dim, 2) *
        gradrsf_QuadPts.block(0, 2 * k, gradrsf_QuadPts.rows(), 2).transpose());
    // Loop over local shape functions/dofs to compute the value of
    // the finite element function in the quadrature point
    Eigen::Matrix<dof_t, Eigen::Dynamic, 1> grad_uh_val(world_dim, 1);
    grad_uh_val.setZero();
    for (int j = 0; j < num_LSF; ++j) {
      grad_uh_val += (dofs[j] * trf_grad.col(j));
    }
    // form the difference
    for (int l = 0; l < world_dim; ++l) {
      grad_uh_val[l] -= vfval[k][l];
    }
    // sum the quared modulus weighted with quadrature weight and metric factor
    sum += qr_Weights[k] * determinants[k] * grad_uh_val.squaredNorm();
    // clang-format on
    SWITCHEDSTATEMENT(ctrl_, kout_comp,
                      std::cout << "\t @ " << (mapped_qpts.col(k)).transpose()
                                << ": grad uh = [" << grad_uh_val.transpose()
                                << "[, vf = [" << vfval[k][0] << ' '
                                << vfval[k][1] << "]" << std::endl);
    // clang-format on
  }
  SWITCHEDSTATEMENT(ctrl_, kout_comp,
                    std::cout << "\t loc norm^2 = " << sum << std::endl);
  return sum;
}

}  // namespace lf::uscalfe

#endif

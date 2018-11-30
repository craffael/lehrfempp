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

namespace lf::fe {
/**
 * @brief Class for the local computation of the L2 norm of the difference of
 *        a Lagrangian finite element function and a generic function
 *
 * @tparam FUNCTOR type for a generic function
 *
 * ### Template parameter type requirements
 *
 * - FUNCTOR must behave like `std::function<SCALAR(Eigen::DenseBase)>` where
 *   SCALAR is a field type.

 *
 * This class complies with the type requirements for the template argument
 * LOCAL_NORM_COMP of the function lf::fe::L2NormDifference. Its main method
 * is the evaluation operator.
 *
 * ### Use case
 *
 * ~~~
 *   LocalL2NormDifference loc_comp(
 *     fe_space.ShapeFunctionLayout(kTria), fe_space.ShapeFunctionLayout(kQuad),
 *      [](auto x) { return (x[0] * (1.0 - x[1])); }, quad_order);
 * ~~~
 */
template <typename FUNCTOR>
class LocalL2NormDifference : public LocCompLagrFEPreprocessor {
 public:
  /** @brief standard constructors */
  /** @{ */
  LocalL2NormDifference(const LocalL2NormDifference &) = delete;
  LocalL2NormDifference(LocalL2NormDifference &&) noexcept = default;
  LocalL2NormDifference &operator=(const LocalL2NormDifference &) = delete;
  LocalL2NormDifference &operator=(LocalL2NormDifference &&) = default;
  /** @} */

  /** @brief Constructor
   *
   * @param fe_tria_p pointer to layout description for reference shape
   * functions on triangular cells
   * @param fe_quad_p pointer to layout description for reference shape
   * functions on quadrilateral cells
   * @param u object providing scalar-valued function
   * @param loc_quad_order desired order of local quadrature, default value = 0.
   *        If = 0, the quadrature order is determined from the polynomial
   *        degree of the reference shape functions.
   *
   * @note NULL pointers amy be passed, if evaluations for the
   *       corresponding cell type are never requested.
   */
  LocalL2NormDifference(
      std::shared_ptr<const ScalarReferenceFiniteElement<double>> fe_tria_p,
      std::shared_ptr<const ScalarReferenceFiniteElement<double>> fe_quad_p,
      FUNCTOR u, lf::quad::quadOrder_t loc_quad_order = 0)
      : LocCompLagrFEPreprocessor(fe_tria_p, fe_quad_p, loc_quad_order),
        u_(u) {}
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

 private:
  /** @brief a handle to the generic function with which to compare the finite
   * element function */
  FUNCTOR u_;

 public:
  /** @brief output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_cell = 8;
  static const unsigned int kout_comp = 16;
};

// Initialization of control variable
template <typename FUNCTOR>
unsigned int LocalL2NormDifference<FUNCTOR>::ctrl_ = 0;

template <typename FUNCTOR>
template <typename DOFVECTOR>
double LocalL2NormDifference<FUNCTOR>::operator()(const lf::mesh::Entity &cell,
                                                  const DOFVECTOR &dofs) {
  using dof_t = typename DOFVECTOR::value_type;
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  LF_ASSERT_MSG(dofs.size() >= num_rsf(ref_el), "Dof vector too short");
  // Obtain cell--type depend invariant data
  const Eigen::MatrixXd &qr_Points{qr_points(ref_el)};
  const Eigen::VectorXd &qr_Weights{qr_weights(ref_el)};
  const size_type qr_NumPts = qr_num_pts(ref_el);
  const size_type num_LSF = num_rsf(ref_el);
  const auto rsf_QuadPts{rsf_at_quadpts(ref_el)};

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
  // Loop over quadrature points
  for (int k = 0; k < qr_NumPts; ++k) {
    // Comparison function value at quadrature point
    const auto uval = u_(mapped_qpts.col(k));
    // Value of finite element function at quadrature point
    // Loop over local shape functions/dofs to compute the value of
    // the finite element function in the quadrature point
    dof_t uh_val = 0.0;
    for (int j = 0; j < num_LSF; ++j) {
      uh_val += rsf_QuadPts(j, k) * dofs[j];
    }
    // form the difference
    uh_val -= uval;
    // sum the quared modulus weighted with quadrature weight and metric factor
    sum += qr_Weights[k] * determinants[k] * std::fabs(uh_val * uh_val);
    SWITCHEDSTATEMENT(ctrl_, kout_comp,
                      std::cout << "\t @ " << (mapped_qpts.col(k)).transpose()
                                << ": uh = " << uh_val << ", u = " << uval
                                << std::endl);
  }
  SWITCHEDSTATEMENT(ctrl_, kout_comp,
                    std::cout << "\t loc norm^2 = " << sum << std::endl);
  return sum;
}

// **********************************************************************

/**
 * @brief Class for the local computation of the L2 norm of the difference of
 *        the _gradient_ of a Lagrangian finite element function and a
 *        generic vector field
 *
 * @tparam VEC_FUNC type for a generic vector field
 *
 * ### Template parameter type requirements
 *
 * - VEC_FUNC must behave like `std::function<VECTOR(Eigen::DenseBase)>` where
 *   VECTOR is a vector type like a vector of Eigen or std::vector.
 *
 * This class complies with the type requirements for the template argument
 * LOCAL_NORM_COMP of the function lf::fe::L2NormDifference. Its main method
 * is the evaluation operator returning the **squared** local norm.
 *
 *
 * ### Use case
 *
 * ~~~
 * todo
 * ~~~
 */
template <typename VEC_FUNC>
class LocL2GradientFEDifference : public LocCompLagrFEPreprocessor {
 public:
  /** @brief type of return value of generic function */
  using vecval_t = typename std::invoke_result<VEC_FUNC, Eigen::Vector2d>::type;

  /** defgroup stdc
   * @brief standard constructors
   * @{ */
  LocL2GradientFEDifference(const LocL2GradientFEDifference &) = delete;
  LocL2GradientFEDifference(LocL2GradientFEDifference &&) noexcept = default;
  LocL2GradientFEDifference &operator=(const LocL2GradientFEDifference &) =
      delete;
  LocL2GradientFEDifference &operator=(LocL2GradientFEDifference &&) = default;
  /** @} */

  /** @brief Constructor
   *
   * @param fe_tria_p reference to layout description for reference shape
   * functions on triangular cells
   * @param fe_quad_p reference to layout description for reference shape
   * functions on quadrilateral cells
   * @param vecfield object providing a generic vector field
   * @param loc_quad_order desired order of local quadrature, default value = 0.
   *        If = 0, the quadrature order is determined from the polynomial
   *        degree of the reference shape functions.
   *
   * @note NULL pointers amy be passed, if evaluations for the
   *       corresponding cell type are never requested.
   */
  LocL2GradientFEDifference(
      std::shared_ptr<const ScalarReferenceFiniteElement<double>> fe_tria_p,
      std::shared_ptr<const ScalarReferenceFiniteElement<double>> fe_quad_p,
      VEC_FUNC vecfield, lf::quad::quadOrder_t loc_quad_order = 0)
      : LocCompLagrFEPreprocessor(fe_tria_p, fe_quad_p, loc_quad_order),
        vecfield_(vecfield) {}
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
   * Actual computation of the element matrix based on numerical quadrature and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying Lagrangian finite element spaces: for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials of degree 2p ("moderate overkill quadrature")
   *
   * ### Template parameter type requirements
   *
   * - DOFVECTOR must be compatible with std::vector of a SCALAR type
   */
  template <typename DOFVECTOR>
  double operator()(const lf::mesh::Entity &cell, const DOFVECTOR &dofs);

 private:
  /** @brief a handle to the generic function with which to compare the finite
   * element function */
  VEC_FUNC vecfield_;

 public:
  /** @brief output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_cell = 8;
  static const unsigned int kout_comp = 16;
};

// Initialization of control variable
template <typename VEC_FUNC>
unsigned int LocL2GradientFEDifference<VEC_FUNC>::ctrl_ = 0;

template <typename VEC_FUNC>
template <typename DOFVECTOR>
double LocL2GradientFEDifference<VEC_FUNC>::operator()(
    const lf::mesh::Entity &cell, const DOFVECTOR &dofs) {
  using dof_t = typename DOFVECTOR::value_type;
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  LF_ASSERT_MSG(dofs.size() >= num_rsf(ref_el), "Dof vector too short");

  // Physical dimension of the cell
  const dim_t world_dim = geo_ptr->DimGlobal();

  // Obtain cell--type depend but cell-independent data
  const Eigen::MatrixXd &qr_Points{qr_points(ref_el)};
  const Eigen::VectorXd &qr_Weights{qr_weights(ref_el)};
  const size_type qr_NumPts = qr_num_pts(ref_el);
  const size_type num_LSF = num_rsf(ref_el);
  const auto gradrsf_QuadPts{grad_rsf_at_quadpts(ref_el)};

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
  // Loop over quadrature points
  for (int k = 0; k < qr_NumPts; ++k) {
    // Value of the vector field at quadrature point
    const auto vfval = vecfield_(mapped_qpts.col(k));
    LF_ASSERT_MSG(vfval.size() == world_dim, "Vector length mismatch");

    // Value of the gradient of the finite element function at quadrature point
    // Transformed gradients
    const auto trf_grad(JinvT.block(0, 2 * k, world_dim, 2) *
                        gradrsf_QuadPts[k]);
    // Loop over local shape functions/dofs to compute the value of
    // the finite element function in the quadrature point
    Eigen::Matrix<dof_t, Eigen::Dynamic, 1> grad_uh_val(world_dim, 1);
    grad_uh_val.setZero();
    for (int j = 0; j < num_LSF; ++j) {
      grad_uh_val += (dofs[j] * trf_grad.col(j));
    }
    // form the difference
    for (int l = 0; l < world_dim; ++l) {
      grad_uh_val[l] -= vfval[l];
    }
    // sum the quared modulus weighted with quadrature weight and metric factor
    sum += qr_Weights[k] * determinants[k] * grad_uh_val.squaredNorm();
    // clang-format on
    SWITCHEDSTATEMENT(ctrl_, kout_comp,
                      std::cout << "\t @ " << (mapped_qpts.col(k)).transpose()
                                << ": grad uh = [" << grad_uh_val.transpose()
                                << "[, vf = [" << vfval[0] << ' ' << vfval[1]
                                << "]" << std::endl);
    // clang-format on
  }
  SWITCHEDSTATEMENT(ctrl_, kout_comp,
                    std::cout << "\t loc norm^2 = " << sum << std::endl);
  return sum;
}

}  // namespace lf::fe

#endif

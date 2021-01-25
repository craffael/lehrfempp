/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

#ifndef LF_LINFE_H
#define LF_LINFE_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Classes supporting the linear finite element discretization of the
 *        Laplacian
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <iostream>

#include "uscalfe.h"

namespace lf::uscalfe {
/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Computing the element matrix for the (negative) Laplacian
 *         and linear finite elements.
 *
 * The main purpose of this class is to compute the element matrix for
 * the Laplacian on affine triangles or bilinearly mapped quadrilaterals.
 * These element matrices are provided by the `Eval()` method.
 *
 * @note the `Eval()` method will always return a _reference_ to a 4x4 matrix
 * also for triangles. In this case the last row and column must be ignored.
 *
 * This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 * @ref AssembleMatrixLocally().
 */
class LinearFELaplaceElementMatrix {
 public:
  using ElemMat = Eigen::Matrix<double, 4, 4>;

  /**
   * @brief Idle constructor
   */
  LinearFELaplaceElementMatrix() = default;
  /**
   * @brief All cells are considered active in the default implementation
   */
  [[nodiscard]] virtual bool isActive(const lf::mesh::Entity & /*cell*/) const {
    return true;
  }
  /**
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a 4x4 matrix, containing the element matrix. The bottom row/column
   *         is not used in the case of a triangle.
   */
  [[nodiscard]] ElemMat Eval(const lf::mesh::Entity &cell) const;

 private:
  /// quadrature points on reference quadrilateral
  const double kSqrt3 = 1.0 / std::sqrt(3.0);
  const double kZeta_0 = 0.5 - 0.5 * kSqrt3;
  const double kZeta_1 = 0.5 + 0.5 * kSqrt3;
  const std::array<Eigen::Vector2d, 4> kQuadPoints{
      Eigen::Vector2d{kZeta_0, kZeta_0}, Eigen::Vector2d{kZeta_0, kZeta_1},
      Eigen::Vector2d{kZeta_1, kZeta_0}, Eigen::Vector2d{kZeta_1, kZeta_1}};
  // Gradients of refrence shape functions in rows of a matrix
  static Eigen::Matrix<double, 4, 2> DervRefShapFncts(
      const Eigen::Vector2d &xh);
  // Barycenter of triangle
  const Eigen::Matrix<double, 2, 1> kTriabc{
      (Eigen::Matrix<double, 2, 1>() << (1.0 / 3), (1.0 / 3)).finished()};
  // Point in reference square
  const Eigen::Matrix<double, 2, 1> kQuadpt{
      (Eigen::Matrix<double, 2, 1>() << 0.7, 0.3).finished()};

 public:
  /**
   * @brief Used by LinearFELaplaceElementMatrix to log additional (debug) info.
   */
  static std::shared_ptr<spdlog::logger> &Logger();
};

/**
 * @ingroup entity_vector_provider
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Class for computation of local load vector for linear finite
 * elements.
 *
 * @tparam FUNCTOR object compatible with a \ref mesh_function type. This means
 *         that is must supply an evaluation operator of signature
 * ~~~
 * std::function<auto(const lf::mesh::Entity &cell,const Eigen::MatrixXd &)>
 * ~~~
 *         which supplies the source function, see @ref
 * mesh::utils::MeshFunctionGlobal.
 *
 * Computations employ edge midpoint quadrature.
 *
 * @note The element vector returned by the `Eval()` method will always
 * have length 4 also for triangles.
 *
 * This class complies with the requirements for the template parameter
 * `ENTITY_VECTOR_PROVIDER` of the function AssembleVectorLocally().
 *
 * TODO: Adjust size of vector returned.
 *
 * #### Logger
 * This class logs additional information to \ref
 * LinearFeLocalLoadVectorLogger().
 * See \ref loggers for more information.
 */
template <typename SCALAR, typename FUNCTOR>
class LinearFELocalLoadVector {
 public:
  using ElemVec = Eigen::Matrix<SCALAR, 4, 1>;

  /** @brief Constructor storing the right hand side function */
  explicit LinearFELocalLoadVector(FUNCTOR f) : f_(f) {}
  /** @brief Default implement: all cells are active */
  [[nodiscard]] virtual bool isActive(const lf::mesh::Entity & /*cell*/) const {
    return true;
  }
  /**
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   *
   * The implementation uses simple edge midpoint based quadrature and an
   * approximation of the volume of a cell just using the integration element at
   * the barycenter.
   */
  [[nodiscard]] ElemVec Eval(const lf::mesh::Entity &cell) const;

 private:
  /** `f_(x)` where `x` is a 2D vector provides the evaluation of the source
   * function */
  FUNCTOR f_;
};

/**
 * @brief Used by LinearFELocalLoadVector to log additional (debug) information
 * during vector assembly.
 */
std::shared_ptr<spdlog::logger> &LinearFeLocalLoadVectorLogger();

// TODO(craffael) remove const once
// https://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <typename SCALAR, typename FUNCTOR>
typename LinearFELocalLoadVector<SCALAR, FUNCTOR>::ElemVec
LinearFELocalLoadVector<SCALAR, FUNCTOR>::Eval(
    const lf::mesh::Entity &cell) const {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  const Eigen::MatrixXd &ref_el_corners(ref_el.NodeCoords());
  // World coordinates of vertices
  // const Eigen::MatrixXd vertices{geo_ptr->Global(ref_el_corners)};
  // Debugging output
  SPDLOG_LOGGER_TRACE(LinearFeLocalLoadVectorLogger(), "{}, shape = \n{}",
                      ref_el, geo_ptr->Global(ref_el_corners));

  // Midpoints of edges in the reference cell
  Eigen::MatrixXd ref_mp(2, 4);
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      // clang-format off
      ref_mp << 0.5,0.5,0.0,-1.0,
                0.0,0.5,0.5,-1.0;
      // clang-format on
      break;
    }
    case lf::base::RefEl::kQuad(): {
      // clang-format off
      ref_mp << 0.5,1,0.5,0.0,
                0.0,0.5,1.0,0.5;
      // clang-format on
      break;
    }
    default: {
      LF_ASSERT_MSG(false, "Illegal entity type!");
      break;
    }
  }  // end switch

  const double area = lf::geometry::Volume(*geo_ptr);

  // Vector for returning element vector
  ElemVec elem_vec = ElemVec::Zero();
  // Run over the midpoints of edges and fetch values of the source function
  // there
  auto fvals = f_(cell, ref_mp);

  for (int k = 0; k < num_nodes; k++) {
    elem_vec[k] += 0.5 * fvals[k];
    elem_vec[(k + 1) % num_nodes] += 0.5 * fvals[k];
  }
  SPDLOG_LOGGER_TRACE(LinearFeLocalLoadVectorLogger(), "element vector = {}",
                      elem_vec.head(num_nodes).transpose());

  return (area / num_nodes) * elem_vec;
}

}  // namespace lf::uscalfe

#endif

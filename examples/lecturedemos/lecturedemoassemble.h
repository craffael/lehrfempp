#ifndef LF_LD_ASS_H
#define LF_LD_ASS_H
/**
 * @file
 * @brief Functions for simple LehrFEM++ demos + sample codes
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include "lecturedemo.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

namespace lecturedemo {

/** @brief Computing the element matrix for the (negative) Laplacian
 *         and linear finite elements on **triangles**
 *
 * The main purpose of this class is to compute the element matrix for
 * the Laplacian on affine triangles. These element matrices are provided
 * by the `Eval()` method.
 *
 * This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 * @ref lf::assemble::AssembleMatrixLocally().
 */
/* SAM_LISTING_BEGIN_1 */
class LinFELaplaceElemMatProvider {
 public:
  /**
   * @brief Idle constructor
   */
  LinFELaplaceElemMatProvider() = default;
  /**
   * @brief All cells are considered active in the default implementation
   */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /**
   * @brief main routine for the computation of element matrices
   *
   * @param tria reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a 3x3 matrix, containing the element matrix.
   */
  Eigen::Matrix<double, 3, 3> Eval(const lf::mesh::Entity &tria);
};
/* SAM_LISTING_END_1 */

/** @brief Computing the edge mass matrix for linear finite elements
 *
 * This class complies with the requirements for the type
 * `ENTITY_MATRIX_PROVIDER` given as a template parameter to define an
 * incarnation of the function
 * @ref lf::assemble::AssembleMatrixLocally().
 */
/* SAM_LISTING_BEGIN_2 */
class LinFEMassEdgeMatProvider {
 public:
  /**
   * @brief Idle constructor
   */
  LinFEMassEdgeMatProvider(lf::mesh::utils::CodimMeshDataSet<bool> &bd_flags)
      : bd_flags_(bd_flags) {}
  /**
   * @brief All cells are considered active in the default implementation
   */
  virtual bool isActive(const lf::mesh::Entity &edge) {
    return bd_flags_(edge);
  }
  /**
   * @brief main routine for the computation of element matrices
   *
   * @param edge reference to an edge of the mesh
   * @return a 2x2 matrix, containing the edge mass matrix.
   */
  Eigen::Matrix2d Eval(const lf::mesh::Entity &edge);

 private:
  lf::mesh::utils::CodimMeshDataSet<bool> &bd_flags_;
};
/* SAM_LISTING_END_2 */

/** @brief Class for computation of local load vector for linear finite
 * elements on triangles.
 *
 * @tparam FUNCTOR object with an evaluation operator of signature
 *         std::function<double(const Eigen::Vector2d &)>, which supplies
 *         the source function
 *
 * Computation is based on vertex based quadrature (2D trapezoidal rule)
 *
 * This class complies with the requirements for the template parameter
 * `ENTITY_VECTOR_PROVIDER` of the function AssembleVectorLocally().
 */
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
class LinFEElemVecProvider {
 public:
  /** @brief Constructor storing the right hand side function */
  explicit LinFEElemVecProvider(FUNCTOR f) : f_(f) {}
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /**
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   *
   * The implementation uses simple vertex based quadrature and an approximation
   * of the volume of a cell just using the integration element at the
   * barycenter.
   */
  Eigen::Vector3d Eval(const lf::mesh::Entity &tria);

 private:
  /** `f_(x)` where `x` is a 2D vector provides the evaluation of the source
   * function */
  FUNCTOR f_;
};
/* SAM_LISTING_END_3 */

// clang-format off
/* SAM_LISTING_BEGIN_4 */
  template <typename FUNCTOR>
  Eigen::Vector3d LinFEElemVecProvider<FUNCTOR>::Eval(
     const lf::mesh::Entity &tria) {
    // Throw error in case no triangular cell
    LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
		  "Unsupported cell type " << tria.RefEl());
    // Obtain vertex coordinates of the triangle in a 2x3 matrix
    const auto corners{lf::geometry::Corners(*(tria.Geometry()))};
    const double area_third = lf::geometry::Volume(*(tria.Geometry())) / 3.0;
    LF_ASSERT_MSG((corners.cols() == 3) && (corners.rows() == 2),
		  "Invalid vertex coordinate " << corners.rows() << "x"
		  << corners.cols() << " matrix");
    return Eigen::Vector3d(area_third * f_(corners.col(0)),
			   area_third * f_(corners.col(1)),
			   area_third * f_(corners.col(2)));
  }
/* SAM_LISTING_END_4 */
// clang-format on

/**
 * @brief Driver routine for demos for LehrFEM++ matrix/vector
 * assembly functions 
 */
void lecturedemoassemble();

  /** 
   * @brief Driver code and demonstration of treatment of Dirichlet boundary conditions
   */
  void lecturedemoDirichlet();

}  // namespace lecturedemo

#endif  // LF_LD_ASS_H

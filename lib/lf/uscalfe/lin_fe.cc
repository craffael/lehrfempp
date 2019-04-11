/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation: element matrix for Laplacian discretized
 *        with linear finite elements
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "lin_fe.h"

namespace lf::uscalfe {
// Implementation for LinearFELaplaceElementMatrix
unsigned int LinearFELaplaceElementMatrix::dbg_ctrl{0};

inline Eigen::Matrix<double, 4, 2>
LinearFELaplaceElementMatrix::DervRefShapFncts(const Eigen::Vector2d &xh) {
  // clang-format off
  return (Eigen::Matrix<double, 4, 2>(4, 2) <<
	  xh[1] - 1, xh[0] - 1,
	  1 - xh[1],  -xh[0]  ,
	  xh[1]    , xh[0]    ,
	  -xh[1]   , 1 - xh[0]
	  ).finished();
  // clang-format on
}

LinearFELaplaceElementMatrix::ElemMat LinearFELaplaceElementMatrix::Eval(
    const lf::mesh::Entity &cell) const {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_geo,
                    std::cout << ref_el << ", shape = \n"
                              << vertices << std::endl);

  // 4x4 dense matrix for returning result
  elem_mat_t elem_mat = elem_mat_t::Zero(4, 4);

  // Computations differ depending on the type of the cell
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      // Triangular cell: straight edges assumed, which means that the triangle
      // is an affine image of the unit triangle and that the gradients of the
      // local shape functions (= barycentric coordinate functions) are
      // constant.
      LF_ASSERT_MSG((vertices.cols() == 3) && (vertices.rows() == 2),
                    "Wrong size of vertex matrix!");
      // Set up an auxiliary 3x3-matrix with a leading column 1 and
      // the vertex coordinates in its right 3x2 block
      Eigen::Matrix<double, 3, 3> X;  // temporary matrix
      X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
      X.block<3, 2>(0, 1) = vertices.transpose();
      // The determinant of the auxliriary matrix also supplies the area
      const double area = 0.5 * std::abs(X.determinant());
      SWITCHEDSTATEMENT(dbg_ctrl, dbg_det,
                        std::cout
                            << "Area: " << area << " <-> "
                            << 0.5 * (geo_ptr->IntegrationElement(kTriabc))
                            << std::endl);

      // Compute the gradients of the barycentric coordinate functions
      // and store them in the columns of a 2x3 matrix grad_bary_coords
      Eigen::Matrix<double, 2, 3> grad_bary_coords{
          X.inverse().block<2, 3>(1, 0)};
      // Since the gradients are constant local integration is easy
      elem_mat.block<3, 3>(0, 0) =
          area * grad_bary_coords.transpose() * grad_bary_coords;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      LF_ASSERT_MSG((vertices.cols() == 4) && (vertices.rows() == 2),
                    "Wrong size of vertex matrix!");
      // Coefficient matrix for bilinear mapping to actual cell
      // clang-format off
      Eigen::Matrix<double, 2, 4> G(vertices *
                                    (Eigen::Matrix<double, 4, 4>() <<
				     1, -1, -1, 1,
				     0, 1 , 0, -1,
				     0, 0, 0, 1,
				     0, 0, 1, -1)
                                        .finished());
      // clang-format on 
      // Coefficients for the determinant as a linear function in
      // the reference coordinates
      // clang-format off
      Eigen::Vector3d detc{
          (Eigen::Vector3d() <<
	   G(0, 1) * G(1, 2) - G(1, 1) * G(0, 2),
           G(0, 1) * G(1, 3) - G(1, 1) * G(0, 3),
           G(1, 2) * G(0, 3) - G(0, 2) * G(1, 3))
              .finished()};
      // clang-format on
      // Determinant at a point in the reference element
      auto detDPhi = [&detc](const Eigen::Vector2d &xh) -> double {
        return std::abs(detc[0] + detc[1] * xh[0] + detc[2] * xh[1]);
      };
      SWITCHEDSTATEMENT(
          dbg_ctrl, dbg_det,
          std::cout << "Determinant: " << detDPhi(kQuadpt) << " <-> "
                    << (geo_ptr->IntegrationElement(kQuadpt)) << std::endl);

      // Transposed adjunct matrix of Jacobian of transformation
      auto DPhiadj =
          [&G](const Eigen::Vector2d &xh) -> Eigen::Matrix<double, 2, 2> {
        // clang-format off
        return ((Eigen::Matrix<double, 2, 2>() <<
		 G(1, 2) + G(1, 3) * xh[0] , -G(1, 1) - G(1, 3) * xh[1],
		 -G(0, 2) - G(0, 3) * xh[0], G(0, 1) + G(0, 3) * xh[1])
                    .finished());
        // clang-format on
      };
      // Output for debugging
      SWITCHEDSTATEMENT(
          dbg_ctrl, dbg_J,
          std::cout << "InverseTransposedJacobian:\n"
                    << (DPhiadj(kQuadpt) / detDPhi(kQuadpt)) << " <-> "
                    << (geo_ptr->JacobianInverseGramian(kQuadpt)) << std::endl);

      // For a quadrilateral we have to use numerical quadrature based on
      // a 2-2 tensor-product Gauss rule.
      // Sum over quadrature points (weight= 0.5)
      for (int q = 0; q < 4; ++q) {
        // Obtain reference coordinates of current quadrature point
        const Eigen::Vector2d qp(kQuadPoints[q]);
        // Gradients of shape functions in quadrature points stored in
        // the columns of a matrix
        Eigen::Matrix<double, 2, 4> gradmat(DPhiadj(qp) *
                                            DervRefShapFncts(qp).transpose());
        // Update of element matrix by contribution from current quadrature
        // point, whose entries are dot products of gradients
        elem_mat += gradmat.transpose() * gradmat / detDPhi(qp);
      }
      // Scaling with quadrature weights
      elem_mat *= 0.25;
      break;
    }
    default: { LF_ASSERT_MSG(false, "Illegal cell type"); }
  }  // end switch

  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_locmat, lf::base::size_type nv = vertices.cols();
      std::cout << "Element matrix\n"
                << elem_mat.block(0, 0, nv, nv) << std::endl;
      std::cout << "Row sums = " << elem_mat.block(0, 0, nv, nv).colwise().sum()
                << ",\n col sums = "
                << elem_mat.block(0, 0, nv, nv).rowwise().sum().transpose()
                << std::endl);

  // Return the element matrix
  return elem_mat;
}

// No implementation for TEMPLATE LinearFELocalLoadVector here

}  // namespace lf::uscalfe

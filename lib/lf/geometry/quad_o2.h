/**
 * @file
 * @brief Declaration of second-order parametric quadrilaterals
 * @author Anian Ruoss
 * @date 2018-12-27 21:17:17
 * @copyright MIT License
 */

#ifndef QUAD_O2_H
#define QUAD_O2_H

#include "geometry_interface.h"

namespace lf::geometry {

/**
 * @brief A second-order quadrilateral in the plane or in 3D space
 *
 * Coordinates \f$ coords = [A, B, C, D, E, F, G, H] \f$ are mapped to the
 * reference element as follows:
 *
 *     (0.0, 1.0) - (0.5, 1.0) - (1.0, 1.0)                  D - G - C
 *          |                         |                      |       |
 *     (0.0, 0.5)                (1.0, 0.5)        ->        H       F
 *          |                         |                      |       |
 *     (0.0, 0.0) - (0.5, 0.0) - (1.0, 0.0)                  A - E - B
 *
 */
class QuadO2 : public Geometry {
 public:
  /**
   * @brief Constructor building quadrilateral from vertex/midpoint coordinates
   * @param coords w x 8 matrix, w = world dimension, whose columns contain the
   *        world coordinates of the vertices/midpoints
   */
  explicit QuadO2(Eigen::Matrix<double, Eigen::Dynamic, 8> coords);

  [[nodiscard]] dim_t DimLocal() const override { return 2; }
  [[nodiscard]] dim_t DimGlobal() const override { return coords_.rows(); }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kQuad();
  }

  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd &local) const override;
  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd &local) const override;
  [[nodiscard]] Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd &local) const override;
  [[nodiscard]] Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd &local) const override;

  /** @copydoc lf::geometry::Geometry::SubGeometry() */
  [[nodiscard]] std::unique_ptr<Geometry> SubGeometry(dim_t codim,
                                                      dim_t i) const override;

  /**
   * @copydoc lf::geometry::Geometry::ChildGeometry()
   *
   * For a detailed description of the indexing of the vertices of child
   * entities see `Refinement.xoj`
   */
  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern &ref_pat, lf::base::dim_t codim) const override;

 private:
  /**
   * @brief Coordinates of the 8 vertices/midpoints, stored in matrix columns
   */
  Eigen::Matrix<double, Eigen::Dynamic, 8> coords_;

  /*
   * QuadO2 is parametrized by:
   *    alpha_ + beta_ * [x1, x2] + gamma_ * [x1^2, x2^2] + delta_ * [x1 * x2]
   *    + epsilon * [x1^2 * x2, x1 * x2^2]
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> beta_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> gamma_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> epsilon_;

  /* Coefficients for efficient evaluation of Jacobian() */
  Eigen::Matrix<double, Eigen::Dynamic, 2> gamma_x_2_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> epsilon_x_2_;
};
}  // namespace lf::geometry

#endif /* QUAD_O2_H */

/**
 * @file
 * @brief Declaration of second-order parametric segments
 * @author Anian Ruoss
 * @date   2018-11-18 19:02:17
 * @copyright MIT License
 */

#ifndef SEGMENT_O2_H
#define SEGMENT_O2_H

#include "geometry_interface.h"

namespace lf::geometry {

/**
 * @brief A second-order segment in the plane or in 3D space.
 *
 * Coordinates \f$ coords = [A, B, C] \f$ are mapped to the reference element as
 * follows:
 *
 *     (0.0) - (0.5) - (1.0)        ->        A - C - B
 *
 */
class SegmentO2 : public Geometry {
 public:
  /**
   * @brief Constructor building segment from vertex/midpoint coordinates
   * @param coords w x 3 matrix, w = world dimension, whose columns contain the
   *        world coordinates of the vertices/midpoints
   */
  explicit SegmentO2(Eigen::Matrix<double, Eigen::Dynamic, 3> coords);

  dim_t DimLocal() const override { return 1; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kSegment(); }

  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;

  /** @copydoc lf::geometry::Geometry::SubGeometry() */
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /**
   * @copydoc lf::geometry::Geometry::ChildGeometry()
   *
   * For a detailed description of the indexing of the vertices of child
   * entities see `Refinement.xoj`
   */
  std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pat, base::dim_t codim) const override;

 private:
  /**
   * @brief Coordinates of the 3 vertices/midpoints, stored in matrix columns
   */
  Eigen::Matrix<double, Eigen::Dynamic, 3> coords_;

  /*
   * SegmentO2 is parametrized by:
   *    alpha_ * x^2 + beta_ * x + gamma_
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> gamma_;

  /*
   * Coefficients for efficient evaluation of JacobianInverseGramian() and
   * IntegrationElement()
   */
  double alpha_squared_;
  double alpha_beta_;
  double beta_squared_;
};

}  // namespace lf::geometry

#endif /* SEGMENT_O2_H */

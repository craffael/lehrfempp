/**
 * @file
 * @brief Declaration of second-order parametric triangles
 * @author Anian Ruoss
 * @date   2018-12-05 22:34:17
 * @copyright MIT License
 */

#ifndef TRIA_O2_H
#define TRIA_O2_H

#include "geometry_interface.h"

namespace lf::geometry {

/**
 * @brief A second-order triangle in the plane or in 3D space.
 *
 * Coordinates \f$ coords = [A, B, C, D, E, F] \f$ are mapped to the reference
 * element as follows:
 *
 *     (0.0, 1.0)                                            C
 *          |     \                                          | \
 *     (0.0, 0.5)   (0.5, 0.5)                     ->        F   E
 *          |                  \                             |     \
 *     (0.0, 0.0) - (0.5, 0.0) - (1.0, 0.0)                  A - D - B
 *
 */
class TriaO2 : public Geometry {
 public:
  explicit TriaO2(Eigen::Matrix<double, Eigen::Dynamic, 6> coords);

  dim_t DimLocal() const override { return 2; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kTria(); }
  Eigen::MatrixXd Global(const Eigen::MatrixXd &local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd &local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd &local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd &local) const override;
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /**
   * @brief creation of child geometries as specified in refinement pattern
   *
   * For a detailed description of the indexing of the vertices of child
   * triangles see `Refinement.xoj`.
   *
   * @sa lf:refinement::Hybrid2DRefinementPattern
   */
  std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern &ref_pat, lf::base::dim_t codim) const override;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 6> coords_;
  // alpha_ + beta_ * [x1, x2] + gamma_ * [x1^2, x2^2] + delta_ * [x1*x2]
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> beta_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> gamma_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> delta_;
  // coefficients for Jacobian
  Eigen::Matrix<double, Eigen::Dynamic, 2> gamma_x_2_;
};
}  // namespace lf::geometry

#endif /* TRIA_O2_H */

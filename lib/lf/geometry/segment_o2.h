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
 * @brief A curved edge parametrized by means of polynomial of degree 2 defined
 * by the location of its two endpoints and its midpoint
 */
class SegmentO2 : public Geometry {
 public:
  explicit SegmentO2(Eigen::Matrix<double, Eigen::Dynamic, 3> coords)
      : coords_(std::move(coords)),
        alpha_(),
        beta_(),
        gamma_(),
        alpha_squared_(),
        alpha_beta_(),
        beta_squared_() {
    const Eigen::VectorXd& vtx0 = coords_.col(0);
    const Eigen::VectorXd& vtx1 = coords_.col(1);
    const Eigen::VectorXd& midp = coords_.col(2);

    // polynomial of degree 2: alpha * x^2 + beta * x + gamma
    alpha_ = 2. * (vtx1 + vtx0) - 4. * midp;
    beta_ = 4. * midp - 3. * vtx0 - vtx1;
    gamma_ = vtx0;

    // coefficients for JacobianInverseGramian and IntegrationElement
    alpha_squared_ = alpha_.squaredNorm();
    alpha_beta_ = alpha_.dot(beta_);
    beta_squared_ = beta_.squaredNorm();
  }

  dim_t DimLocal() const override { return 1; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kSegment(); }
  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /** @brief creation of child geometry for the sake of mesh refinement
   *
   * @see Geometry::ChildGeometry()
   * @see RefinementPattern
   * @param ref_pat three refinement patterns are supported
   * - rp_nil: empty refinement
   * - rp_copy: just copies the geometry information of the segment
   * - rp_split: split edge in the middle.
   * @param codim _relative_ codimension of children whose shape is requested
   */
  std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pat, base::dim_t codim) const override;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 3> coords_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> gamma_;
  double alpha_squared_;
  double alpha_beta_;
  double beta_squared_;
};

}  // namespace lf::geometry

#endif /* SEGMENT_O2_H */

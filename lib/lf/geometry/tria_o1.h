#ifndef __b84bb391fa744cdd9159293bfcb5b311
#define __b84bb391fa744cdd9159293bfcb5b311

#include "geometry_interface.h"

namespace lf::geometry {

/** @brief Asserting non-degenerate shape of a flat triangle
 *
 * @param coords w x 3 Eigen matrix whose columns contain the vertex coordinates
 *        of the quadrilateral, w = world dimension
 * @param tol relative tolerance for numerical tests of equality with zero
 *
 * Terminates execution in case degenerate shape is detected.
 */
bool assertNonDegenerateTriangle(
    const Eigen::Matrix<double, Eigen::Dynamic, 3>& coords,
    double tol = 1.0E-8);

/**
 * @brief An affine triangle in the plane or in 3D space
 */
class TriaO1 : public Geometry {
 public:
  explicit TriaO1(Eigen::Matrix<double, Eigen::Dynamic, 3> coords);

  [[nodiscard]] dim_t DimLocal() const override { return 2; }
  [[nodiscard]] dim_t DimGlobal() const override { return coords_.rows(); }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override {
    return jacobian_.replicate(1, local.cols());
  }
  [[nodiscard]] Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override {
    return jacobian_inverse_gramian_.replicate(1, local.cols());
  }
  [[nodiscard]] Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override {
    return Eigen::VectorXd::Constant(local.cols(), integrationElement_);
  }
  [[nodiscard]] std::unique_ptr<Geometry> SubGeometry(dim_t codim,
                                                      dim_t i) const override;

  /**
   * @brief creation of child geometries as specified in refinement pattern
   *
   * For a detailed description of the indexing of the vertices of child
   * triangles see `Refinement.xoj`.
   *
   * @sa lf:refinement::Hybrid2DRefinementPattern
   */
  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pat, lf::base::dim_t codim) const override;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 3> coords_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> jacobian_;
  Eigen::Matrix<double, Eigen::Dynamic, 2> jacobian_inverse_gramian_;
  double integrationElement_;
};
}  // namespace lf::geometry

#endif  // __b84bb391fa744cdd9159293bfcb5b311

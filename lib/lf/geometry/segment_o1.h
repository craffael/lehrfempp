#ifndef __0ec30f6fc5554ec9bea4d9b83e83ea42
#define __0ec30f6fc5554ec9bea4d9b83e83ea42

#include "geometry_interface.h"

namespace lf::geometry {

  /**
   * @brief A straight edge defined by the location of its two endpoints
   */
class SegmentO1 : public Geometry {
 public:
  explicit SegmentO1(Eigen::Matrix<double, Eigen::Dynamic, 2> coords)
      : coords_(std::move(coords)) {}

  dim_t DimLocal() const override { return 1; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kSegment(); }
  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const ::Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /** @brief creation of child geometry for the sake of mesh refinement 
   *
   * @see Geometry::ChildGeometry()
   * @see RefinementPattern
   * @param ref_pattern three refinement patterns are supported
   * - rp_copy: just copies the geometry information of the segment
   * - rp_split, rp_regular: split edge in the middle.
   * @param anchor not used
   * @param selector = 0 for splitting refinement pattern returns half edge
   * adjacent to endpoint 0, =1 returns the other half.
   */
  virtual std::vector<std::unique_ptr<Geometry>>
  ChildGeometry(int ref_pattern,int anchor, int selector) const override;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 2> coords_;
};

}  // namespace lf::geometry

#endif  // __0ec30f6fc5554ec9bea4d9b83e83ea42

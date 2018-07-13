#ifndef __b193e889f4e245959d124b0ded8a3b68
#define __b193e889f4e245959d124b0ded8a3b68

#include "geometry_interface.h"

namespace lf::geometry {

/**
 * @brief Bilinear quadrilateral element shape
 *
 * This element shape arises from subjecting the unit square to
 * a componentwise bilinear mapping. This results in a _non-flat_ quadrilateral
 * with _straight_ edges, whose shape is determined by the location of the
 * four vertices.
 */
class QuadO1 : public Geometry {
 public:
  explicit QuadO1(Eigen::Matrix<double, Eigen::Dynamic, 4> coords);

  dim_t DimLocal() const override { return 2; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kQuad(); }

  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const ::Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;

  /** @copydoc Geometry::SubGeometry() */
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /** @copydoc Geometry::isAffine() */
  bool isAffine() const override { return false; }

  /**
   * @brief creation of child geometries as specified in refinement pattern
   *
   * For a detailed description of the indexing of the vertices of child
   * entities see `Refinement.xoj`.
   */
  std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pat) const override;

 private:
  /** @brief Coordinates of the a four vertices */
  Eigen::Matrix<double, Eigen::Dynamic, 4> coords_;
};

}  // namespace lf::geometry

#endif  // __b193e889f4e245959d124b0ded8a3b68

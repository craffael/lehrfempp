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
 * @brief A second-order triangle in the plane or in 3D space
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
};
}  // namespace lf::geometry

#endif /* TRIA_O2_H */

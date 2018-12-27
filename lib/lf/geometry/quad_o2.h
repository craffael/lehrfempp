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
class QuadO2 : public Geometry {
 public:
  /**
   * @brief Constructor build quadrilateral from vertex coordinates
   * @param coords w x 8 matrix, w = world dimension, whose columns contain the
   *        world coordinates of the vertices
   */
  explicit QuadO2(Eigen::Matrix<double, Eigen::Dynamic, 8> coords);

  dim_t DimLocal() const override { return 2; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kQuad(); }

  Eigen::MatrixXd Global(const Eigen::MatrixXd &local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd &local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd &local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd &local) const override;

  /** @copydoc Geometry::SubGeometry() */
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /**
   * @copydoc lf::geometry::Geometry::ChildGeometry()
   *
   * For a detailed description of the indexing of the vertices of child
   * entities see `Refinement.xoj`
   */
  std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern &ref_pat, lf::base::dim_t codim) const override;

 private:
  /**
   * @brief Coordinates of the 8 vertices/midpoints, stored in matrix columns
   */
  Eigen::Matrix<double, Eigen::Dynamic, 8> coords_;
};
}  // namespace lf::geometry

#endif /* QUAD_O2_H */

/**
 * @file
 * @brief Defines the OcctBrepModel class
 * @author Raffael Casagrande
 * @date   2020-11-06 05:12:41
 * @copyright MIT License
 */

#ifndef __a7dfc983b4ff4ca19ec578d1a86003c5
#define __a7dfc983b4ff4ca19ec578d1a86003c5

#include <Bnd_OBB.hxx>
#include <Geom_Curve.hxx>
#include <TopoDS_Edge.hxx>
#include "lf/brep/interface/interface.h"

namespace lf::brep::occt {

class OcctCurveGeometry final : public interface::BrepGeometry {
 public:
  explicit OcctCurveGeometry(const Bnd_OBB &obb, const TopoDS_Edge &edge);

  OcctCurveGeometry(const OcctCurveGeometry &) = default;
  OcctCurveGeometry(OcctCurveGeometry &&) = default;
  OcctCurveGeometry &operator=(const OcctCurveGeometry &) = default;
  OcctCurveGeometry &operator=(OcctCurveGeometry &&) = default;
  virtual ~OcctCurveGeometry() = default;

  [[nodiscard]] base::dim_t DimGlobal() const override { return 3; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd &local) const override;
  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd &local) const override;
  [[nodiscard]] std::pair<Eigen::VectorXd, Eigen::MatrixXd> Project(
      const Eigen::MatrixXd &global) const override;
  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd &global) const override;
  [[nodiscard]] std::vector<bool> IsInside(
      const Eigen::MatrixXd &local) const override;

 private:
  TopoDS_Edge edge_;
  Bnd_OBB obb_;
  opencascade::handle<Geom_Curve> curve_;
  double umin_;
  double umax_;
};

}  // namespace lf::brep::occt

#endif  // __a7dfc983b4ff4ca19ec578d1a86003c5

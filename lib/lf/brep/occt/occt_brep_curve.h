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

class OcctBrepCurve final : public interface::BrepGeometry {
 public:
  explicit OcctBrepCurve(TopoDS_Edge &&edge);

  OcctBrepCurve(const OcctBrepCurve &) = default;
  OcctBrepCurve(OcctBrepCurve &&) = default;
  OcctBrepCurve &operator=(const OcctBrepCurve &) = default;
  OcctBrepCurve &operator=(OcctBrepCurve &&) = default;
  virtual ~OcctBrepCurve() = default;

  [[nodiscard]] base::dim_t DimGlobal() const override { return 3; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd &local) const override;

  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd &local) const override;

  [[nodiscard]] std::pair<double, Eigen::VectorXd> Project(
      const Eigen::VectorXd &global) const override;
  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd &global) const override;

  [[nodiscard]] bool IsInside(const Eigen::VectorXd &local) const override;

  [[nodiscard]] Eigen::VectorXd Periods() const override;

  // OCCT specific member functions:

  [[nodiscard]] const TopoDS_Edge &Edge() const { return edge_; }

 private:
  TopoDS_Edge edge_;
  Bnd_OBB obb_;
  opencascade::handle<Geom_Curve> curve_;
  double umin_;
  double umax_;
};

}  // namespace lf::brep::occt

#endif  // __a7dfc983b4ff4ca19ec578d1a86003c5

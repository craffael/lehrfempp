/**
 * @file
 * @brief Implementation of OcctBrepGeometry
 * @author Raffael Casagrande
 * @date   2020-11-06 05:17:46
 * @copyright MIT License
 */

#include "occt_curve_geometry.h"
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

#include "occt_utils.h"

namespace lf::brep::occt {
namespace /* anonymous */ {}  // namespace

OcctCurveGeometry::OcctCurveGeometry(const Bnd_OBB &obb,

                                     const TopoDS_Edge &edge)
    : edge_(edge), obb_(obb) {
  curve_ = BRep_Tool::Curve(edge_, umin_, umax_);
}

Eigen::MatrixXd OcctCurveGeometry::Global(const Eigen::MatrixXd &local) const {
  LF_ASSERT_MSG(local.rows() == 1,
                "local must have exactly one row because this is a curve (1D).")
  Eigen::MatrixXd result(3, local.cols());
  gp_Pnt p;
  for (int i = 0; i < local.cols(); ++i) {
    curve_->D0(local(0, i), p);
    result.col(i) = detail::ToVector(p);
  }
  return result;
}

Eigen::MatrixXd OcctCurveGeometry::Jacobian(
    const Eigen::MatrixXd &local) const {
  LF_ASSERT_MSG(local.rows() == 1,
                "local must have exactly one row because this is a curve (1D).")
  Eigen::MatrixXd result(3, local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    gp_Pnt p;
    gp_Vec v;
    curve_->D1(local(0, i), p, v);
    result.col(i) = detail::ToVector(v);
  }

  return result;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> OcctCurveGeometry::Project(
    const Eigen::MatrixXd &global) const {
  LF_ASSERT_MSG(global.rows() == 2 || global.rows() == 3,
                "global must have 2 or 3 rows.")
  Eigen::VectorXd distance(1, global.cols());
  Eigen::MatrixXd parameters(1, global.cols());

  for (int i = 0; i < global.cols(); ++i) {
    GeomAPI_ProjectPointOnCurve proj(detail::ToPoint(global.col(i)), curve_,
                                     umin_, umax_);
    parameters(0, i) = proj.LowerDistanceParameter();
    distance[i] = proj.LowerDistance();
  }
  return {std::move(distance), std::move(parameters)};
}

std::vector<bool> OcctCurveGeometry::IsInBoundingBox(
    const Eigen::MatrixXd &global) const {
  LF_ASSERT_MSG(global.rows() == 2 || global.rows() == 3,
                "global must have 2 or 3 rows.")

  std::vector<bool> result(global.cols());
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = !obb_.IsOut(detail::ToPoint(global.col(i)));
  }
  return result;
}

std::vector<bool> OcctCurveGeometry::IsInside(
    const Eigen::MatrixXd &local) const {
  LF_ASSERT_MSG(local.rows() == 1,
                "local must have exactly one row because this is a curve (1D).")
  std::vector<bool> result(local.cols());

  for (int i = 0; i < local.cols(); ++i) {
    result[i] = (umin_ - local(0, i) < Precision::PApproximation()) &&
                (local(0, i) - umax_ < Precision::PApproximation());
  }
  return result;
}

}  // namespace lf::brep::occt

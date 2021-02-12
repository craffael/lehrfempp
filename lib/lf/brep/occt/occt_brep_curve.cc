/**
 * @file
 * @brief Implementation of OcctBrepGeometry
 * @author Raffael Casagrande
 * @date   2020-11-06 05:17:46
 * @copyright MIT License
 */

#include "occt_brep_curve.h"

#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

#include "lf/mesh/utils/all_codim_mesh_data_set.h"
#include "occt_details.h"

namespace lf::brep::occt {
namespace /* anonymous */ {}  // namespace

OcctBrepCurve::OcctBrepCurve(TopoDS_Edge &&edge) : edge_(edge) {
  BRepBndLib::AddOBB(edge_, obb_);
  curve_ = BRep_Tool::Curve(edge_, umin_, umax_);

  LF_VERIFY_MSG(!curve_.IsNull(),
                "Unexpected: there is no curve associated with this edge.");
}

Eigen::MatrixXd OcctBrepCurve::Global(const Eigen::MatrixXd &local) const {
  LF_ASSERT_MSG(local.rows() == 1,
                "local must have exactly one row because this is a curve (1D).")
  Eigen::MatrixXd result(3, local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    gp_Pnt p;
    curve_->D0(local(0, i), p);
    result.col(i) = detail::ToVector(p);
  }
  return result;
}

Eigen::MatrixXd OcctBrepCurve::Jacobian(const Eigen::MatrixXd &local) const {
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

std::pair<double, Eigen::VectorXd> OcctBrepCurve::Project(
    const Eigen::VectorXd &global) const {
  LF_ASSERT_MSG(global.rows() == 3, "Unexpected #rows in global.");
  GeomAPI_ProjectPointOnCurve proj(detail::ToPoint(global), curve_);
  Eigen::VectorXd param(1);
  param(0) = proj.LowerDistanceParameter();
  return {proj.LowerDistance(), std::move(param)};
}

std::vector<bool> OcctBrepCurve::IsInBoundingBox(
    const Eigen::MatrixXd &global) const {
  LF_ASSERT_MSG(global.rows() == 3, "global must have 3 rows.")

  std::vector<bool> result(global.cols());
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = !obb_.IsOut(detail::ToPoint(global.col(i)));
  }
  return result;
}

bool OcctBrepCurve::IsInside(const Eigen::VectorXd &local) const {
  LF_ASSERT_MSG(local.rows() == 1, "local should have exactly 1 row.");
  return (umin_ - local(0) < Precision::PApproximation()) &&
         (local(0) - umax_ < Precision::PApproximation());
}

Eigen::VectorXd OcctBrepCurve::Periods() const {
  Eigen::VectorXd result(1);
  if (curve_->IsPeriodic()) {
    result(0) = curve_->Period();
  } else {
    result(0) = 0;
  }
  return result;
}

}  // namespace lf::brep::occt

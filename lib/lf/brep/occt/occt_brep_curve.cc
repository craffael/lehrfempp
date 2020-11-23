/**
 * @file
 * @brief Implementation of OcctBrepGeometry
 * @author Raffael Casagrande
 * @date   2020-11-06 05:17:46
 * @copyright MIT License
 */

#include "occt_brep_curve.h"
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

#include "lf/mesh/utils/all_codim_mesh_data_set.h"
#include "occt_details.h"

#include <BRepBndLib.hxx>

namespace lf::brep::occt {
namespace /* anonymous */ {}  // namespace

OcctBrepCurve::OcctBrepCurve(TopoDS_Edge &&edge) : edge_(edge) {
  BRepBndLib::AddOBB(edge_, obb_);
  curve_ = BRep_Tool::Curve(edge_, umin_, umax_);

  LF_VERIFY_MSG(!curve_.IsNull(),
                "Unexpected: there is no curve associated with this edge.");
}

Eigen::MatrixXd OcctBrepCurve::GlobalMulti(const Eigen::MatrixXd &local) const {
  LF_ASSERT_MSG(local.rows() == 1,
                "local must have exactly one row because this is a curve (1D).")
  Eigen::MatrixXd result(3, local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    result.col(i) = GlobalSingle(local(0, i));
  }
  return result;
}

Eigen::Vector3d OcctBrepCurve::GlobalSingle(double local) const {
  gp_Pnt p;
  curve_->D0(local, p);
  return detail::ToVector(p);
}

Eigen::MatrixXd OcctBrepCurve::JacobianMulti(
    const Eigen::MatrixXd &local) const {
  LF_ASSERT_MSG(local.rows() == 1,
                "local must have exactly one row because this is a curve (1D).")
  Eigen::MatrixXd result(3, local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    result.col(i) = JacobianSingle(local(0, i));
  }
  return result;
}

Eigen::Vector3d OcctBrepCurve::JacobianSingle(double local) const {
  gp_Pnt p;
  gp_Vec v;
  curve_->D1(local, p, v);
  return detail::ToVector(v);
}

std::pair<double, double> OcctBrepCurve::Project(
    const Eigen::Vector3d &global) const {
  GeomAPI_ProjectPointOnCurve proj(detail::ToPoint(global), curve_);
  return {proj.LowerDistance(), proj.LowerDistanceParameter()};
}

std::vector<bool> OcctBrepCurve::IsInBoundingBoxMulti(
    const Eigen::MatrixXd &global) const {
  LF_ASSERT_MSG(global.rows() == 3, "global must have 3 rows.")

  std::vector<bool> result(global.cols());
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = !obb_.IsOut(detail::ToPoint(global.col(i)));
  }
  return result;
}

bool OcctBrepCurve::IsInBoundingBoxSingle(
    const Eigen::Vector3d &global) const {
  return !obb_.IsOut(detail::ToPoint(global));
}

bool OcctBrepCurve::IsInside(double local) const {
  return (umin_ - local < Precision::PApproximation()) &&
         (local - umax_ < Precision::PApproximation());
}

}  // namespace lf::brep::occt

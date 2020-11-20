/**
 * @file
 * @brief Implementation of OcctSurfaceGeometry
 * @author Raffael Casagrande
 * @date   2020-11-20 03:35:04
 * @copyright MIT License
 */

#include "occt_surface_geometry.h"

#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <IntTools_Tools.hxx>

#include "occt_utils.h"

namespace lf::brep::occt {
OcctSurfaceGeometry::OcctSurfaceGeometry(TopoDS_Face&& face)
    : face_(std::move(face)) {
  BRepBndLib::AddOBB(face_, obb_);
  surface_ = BRep_Tool::Surface(face_);
}

Eigen::MatrixXd OcctSurfaceGeometry::Global(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(
      local.rows() == 2,
      "local must have exactly two rows because this is a curve (2d)");
  Eigen::MatrixXd result(3, local.cols());
  gp_Pnt p;
  for (int i = 0; i < local.cols(); ++i) {
    surface_->D0(local(0, i), local(1, i), p);
    result.col(i) = detail::ToVector(p);
  }
  return result;
}

Eigen::MatrixXd OcctSurfaceGeometry::Jacobian(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(false, "not implemented");
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> OcctSurfaceGeometry::Project(
    const Eigen::MatrixXd& global) const {
  LF_ASSERT_MSG(global.rows() == 3, "global must have 3 rows.");
  Eigen::VectorXd distance(global.cols());
  Eigen::MatrixXd parameters(2, global.cols());

  for (int i = 0; i < global.cols(); ++i) {
    GeomAPI_ProjectPointOnSurf proj(detail::ToPoint(global.col(i)), surface_);
    distance[i] = proj.LowerDistance();
    proj.LowerDistanceParameters(parameters(0, i), parameters(1, i));
  }

  return {std::move(distance), std::move(parameters)};
}

std::vector<bool> OcctSurfaceGeometry::IsInBoundingBox(
    const Eigen::MatrixXd& global) const {
  LF_ASSERT_MSG(global.rows() == 3, "global must have 3 rows.");

  std::vector<bool> result(global.cols());
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = !obb_.IsOut(detail::ToPoint(global.col(i)));
  }
  return result;
}

std::vector<bool> OcctSurfaceGeometry::IsInside(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 2,
                "local must have 2 rows because this is a surface (2D).");
  std::vector<bool> result(local.cols());

  for (int i = 0; i < local.cols(); ++i) {
    auto state = IntTools_Tools::ClassifyPointByFace(
        face_, detail::ToPoint2d(local.col(i)));
    result[i] = (state == TopAbs_IN || state == TopAbs_ON);
  }
  return result;
}
}  // namespace lf::brep::occt

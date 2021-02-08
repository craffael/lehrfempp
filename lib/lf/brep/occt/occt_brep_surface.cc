/**
 * @file
 * @brief Implementation of OcctSurfaceGeometry
 * @author Raffael Casagrande
 * @date   2020-11-20 03:35:04
 * @copyright MIT License
 */

#include "occt_brep_surface.h"

#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <IntTools_Tools.hxx>

#include "occt_brep_surface.h"
#include "occt_details.h"

namespace lf::brep::occt {
OcctBrepSurface::OcctBrepSurface(TopoDS_Face&& face) : face_(std::move(face)) {
  BRepBndLib::AddOBB(face_, obb_);
  surface_ = BRep_Tool::Surface(face_);
  LF_VERIFY_MSG(!surface_.IsNull(),
                "Unexpected: there is no surface associated with this face.");
}

Eigen::MatrixXd OcctBrepSurface::GlobalMulti(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(
      local.rows() == 2,
      "local must have exactly two rows because this is a curve (2d)");
  Eigen::MatrixXd result(3, local.cols());

  for (int i = 0; i < local.cols(); ++i) {
    result.col(i) = GlobalSingle(local.col(i));
  }
  return result;
}

Eigen::Vector3d OcctBrepSurface::GlobalSingle(
    const Eigen::Vector2d& local) const {
  return detail::ToVector(surface_->Value(local.x(), local.y()));
}

Eigen::MatrixXd OcctBrepSurface::JacobianMulti(
    const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(3, 2 * local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    result.block<3, 2>(0, 2 * i) = JacobianSingle(local.col(i));
  }
  return result;
}

Eigen::Matrix<double, 3, 2> OcctBrepSurface::JacobianSingle(
    const Eigen::Vector2d& local) const {
  Eigen::Matrix<double, 3, 2> result;
  gp_Pnt p;
  surface_->D1(local(0), local(1), p, *reinterpret_cast<gp_Vec*>(result.data()),
               *reinterpret_cast<gp_Vec*>(result.data() + 3));
  return result;
}

std::pair<double, Eigen::Vector2d> OcctBrepSurface::Project(
    const Eigen::Vector3d& global) const {
  GeomAPI_ProjectPointOnSurf proj(detail::ToPoint(global), surface_);

  Eigen::Vector2d v;
  proj.LowerDistanceParameters(v.x(), v.y());
  return {proj.LowerDistance(), v};
}

std::vector<bool> OcctBrepSurface::IsInBoundingBoxMulti(
    const Eigen::MatrixXd& global) const {
  LF_ASSERT_MSG(global.rows() == 3, "global must have 3 rows.");

  std::vector<bool> result(global.cols());
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = !obb_.IsOut(detail::ToPoint(global.col(i)));
  }
  return result;
}

bool OcctBrepSurface::IsInBoundingBox(const Eigen::Vector3d& global) const {
  return !obb_.IsOut(detail::ToPoint(global));
}

bool OcctBrepSurface::IsInside(const Eigen::Vector2d& local) const {
  auto state =
      IntTools_Tools::ClassifyPointByFace(face_, detail::ToPoint2d(local));
  return (state == TopAbs_IN || state == TopAbs_ON);
}
}  // namespace lf::brep::occt

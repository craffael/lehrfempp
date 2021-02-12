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

Eigen::MatrixXd OcctBrepSurface::Global(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(
      local.rows() == 2,
      "local must have exactly two rows because this is a curve (2d)");
  Eigen::MatrixXd result(3, local.cols());

  for (int i = 0; i < local.cols(); ++i) {
    result.col(i) = detail::ToVector(surface_->Value(local(0, i), local(1, i)));
  }
  return result;
}

Eigen::MatrixXd OcctBrepSurface::Jacobian(
    const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(3, 2 * local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    gp_Pnt p;
    surface_->D1(local(0, i), local(1, i), p,
                 *reinterpret_cast<gp_Vec*>(&result(0, 2 * i)),
                 *reinterpret_cast<gp_Vec*>(&result(0, 2 * i + 1)));
  }
  return result;
}

std::pair<double, Eigen::VectorXd> OcctBrepSurface::Project(
    const Eigen::VectorXd& global) const {
  LF_ASSERT_MSG(global.rows() == 3, "Unexpected number of global points.");
  GeomAPI_ProjectPointOnSurf proj(detail::ToPoint(global), surface_);

  Eigen::Vector2d v;
  proj.LowerDistanceParameters(v.x(), v.y());
  return {proj.LowerDistance(), v};
}

std::vector<bool> OcctBrepSurface::IsInBoundingBox(
    const Eigen::MatrixXd& global) const {
  LF_ASSERT_MSG(global.rows() == 3, "global must have 3 rows.");

  std::vector<bool> result(global.cols());
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = !obb_.IsOut(detail::ToPoint(global.col(i)));
  }
  return result;
}

bool OcctBrepSurface::IsInside(const Eigen::VectorXd& local) const {
  LF_ASSERT_MSG(local.rows() == 2, "Illegal number of rows in local.");
  auto state =
      IntTools_Tools::ClassifyPointByFace(face_, detail::ToPoint2d(local));
  return (state == TopAbs_IN || state == TopAbs_ON);
}
}  // namespace lf::brep::occt

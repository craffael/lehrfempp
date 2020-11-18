/**
 * @file
 * @brief Implementation of OcctBrepModel
 * @author Raffael Casagrande
 * @date   2020-11-06 05:19:08
 * @copyright MIT License
 */

#include "occt_brep_model.h"
#include <BRepBndLib.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>

#include "occt_curve_geometry.h"
#include "occt_utils.h"

namespace lf::brep::occt {
OcctBrepModel::OcctBrepModel(std::string_view filename) {
  BRep_Builder builder;
  if (!BRepTools::Read(shape_, filename.data(), builder)) {
    LF_VERIFY_MSG(false, "Could not open file " << filename);
  }
  LF_VERIFY_MSG(shape_.ShapeType() == TopAbs_COMPOUND,
                "expected a compound type");

  // retrieve all edges with their bounding boxes and save them:

  auto explorer = TopExp_Explorer(shape_, TopAbs_EDGE);
  while (explorer.More()) {
    auto edge = TopoDS::Edge(explorer.Current());

    // check if a same edge (same TShape, same location) already exists, and if
    // this is the case skip it:
    if (std::any_of(edges_.begin(), edges_.end(),
                    [&](const auto& e) { return e.second.IsSame(edge); })) {
      explorer.Next();
      continue;
    }

    Bnd_OBB obb;
    BRepBndLib::AddOBB(edge, obb);
    edges_.emplace_back(std::move(obb), std::move(edge));
    explorer.Next();
  }

  // retrieve all faces with their bounding boxes and save them:
  explorer.Init(shape_, TopAbs_FACE);
  while (explorer.More()) {
    auto face = TopoDS::Face(explorer.Current());
    Bnd_OBB obb;
    BRepBndLib::AddOBB(face, obb);
    faces_.emplace_back(std::move(obb), std::move(face));
    explorer.Next();
  }
}

inline std::pair<std::unique_ptr<interface::BrepGeometry>, Eigen::MatrixXd>
OcctBrepModel::FindGeometry(base::dim_t dim,
                            const Eigen::MatrixXd& global) const {
  LF_VERIFY_MSG(dim > 0 && dim < 3, "dim must be either 1 or 2.");
  LF_ASSERT_MSG(global.cols() > 0, "global must contain at least one column.");
  LF_ASSERT_MSG(global.rows() > 1 && global.rows() <= 3,
                "global must contain 2 or 3 rows.");

  std::vector<gp_Pnt> points;
  points.reserve(global.cols());
  for (int c = 0; c < global.cols(); ++c) {
    points.emplace_back(detail::ToPoint(global.col(c)));
  }

  if (dim == 1) {
    Eigen::MatrixXd result_coord(1, global.cols());

    for (const auto& e : edges_) {
      // 1) Check if all points are inside the bounding box:
      if (std::any_of(points.begin(), points.end(),
                      [&](const auto& p) { return e.first.IsOut(p); })) {
        continue;
      }

      // 2) project point on curve and check if it is within the bounds:
      double start, end;
      auto curve = BRep_Tool::Curve(e.second, start, end);
      LF_VERIFY_MSG(!curve.IsNull(),
                    "Unexpected: there is no curve associated with this edge.");
      for (unsigned i = 0; i < points.size(); ++i) {
        GeomAPI_ProjectPointOnCurve proj(points[i], curve, start, end);

        if (proj.NbPoints() == 0 || proj.LowerDistance() > 1e-5) {
          break;
        }
        result_coord(0, i) = proj.LowerDistanceParameter();
        if (i + 1 == points.size()) {
          return {std::make_unique<OcctCurveGeometry>(e.first, e.second),
                  std::move(result_coord)};
        }
      }
    }

    // if not found, return nullptr:
    return {nullptr, Eigen::MatrixXd()};
  } else {
    LF_ASSERT_MSG(false, "not yet implemented.");
  }
}

}  // namespace lf::brep::occt

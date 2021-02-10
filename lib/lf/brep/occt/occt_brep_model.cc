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
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>

#include "occt_brep_curve.h"
#include "occt_details.h"

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

    if (BRep_Tool::Degenerated(edge)) {
      // This is a degenerated edge (e.g. the pole edge of a sphere) and doesn't
      // have a geometry -> ignore it.
      explorer.Next();
      continue;
    }

    // check if a same edge (same TShape, same location) already exists, and if
    // this is the case skip it:
    if (std::any_of(edges_.begin(), edges_.end(),
                    [&](const auto& e) { return edge.IsSame(e->Edge()); })) {
      explorer.Next();
      continue;
    }

    edges_.push_back(std::make_shared<OcctBrepCurve>(std::move(edge)));
    explorer.Next();
  }

  // retrieve all faces with their bounding boxes and save them:
  explorer.Init(shape_, TopAbs_FACE);
  while (explorer.More()) {
    auto face = TopoDS::Face(explorer.Current());

    if (std::any_of(faces_.begin(), faces_.end(),
                    [&](const auto& f) { return face.IsSame(f->Face()); })) {
      explorer.Next();
      continue;
    }

    faces_.push_back(std::make_shared<OcctBrepSurface>(std::move(face)));
    explorer.Next();
  }
}

std::vector<std::pair<std::shared_ptr<const interface::BrepCurve>, double>>
OcctBrepModel::FindCurves(const Eigen::Vector3d& global) const {
  std::vector<std::pair<std::shared_ptr<interface::BrepCurve const>, double>>
      result;
  for (const auto& e : edges_) {
    if (!e->IsInBoundingBoxSingle(global)) {
      continue;
    }
    auto [dist, param] = e->Project(global);
    if (dist <= 1e-5 && e->IsInside(param)) {
      result.emplace_back(e, param);
    }
  }
  return result;
}

inline std::vector<
    std::pair<std::shared_ptr<const interface::BrepCurve>, Eigen::RowVectorXd>>
OcctBrepModel::FindCurvesMulti(const Eigen::Matrix3Xd& global) const {
  LF_ASSERT_MSG(global.cols() > 0, "global must contain at least one column.");

  std::vector<std::pair<std::shared_ptr<const interface::BrepCurve>,
                        Eigen::RowVectorXd>>
      result;

  Eigen::RowVectorXd result_coord;

  for (const auto& e : edges_) {
    // 1) Check if all points are inside the bounding box:
    if ([&]() {
          for (int i = 0; i < global.cols(); ++i) {
            if (!e->IsInBoundingBoxSingle(global.col(i))) {
              return true;
            }
          }
          return false;
        }()) {
      continue;
      ;
    }

    // 2) project point on curve and check if it is within the bounds:
    for (unsigned i = 0; i < global.cols(); ++i) {
      auto [dist, param] = e->Project(global.col(i));
      if (dist > 1e-5 || !e->IsInside(param)) {
        break;
      }

      if (result_coord.cols() == 0) {
        result_coord.resize(1, global.cols());
      }

      result_coord(0, i) = param;
      if (i + 1 == global.cols()) {
        result.emplace_back(e, std::move(result_coord));
        result_coord = Eigen::RowVectorXd();
      }
    }
  }

  return result;
}

std::vector<
    std::pair<std::shared_ptr<const interface::BrepSurface>, Eigen::Vector2d>>
OcctBrepModel::FindSurfaces(const Eigen::Vector3d& global) const {
  std::vector<
      std::pair<std::shared_ptr<const interface::BrepSurface>, Eigen::Vector2d>>
      result;
  for (const auto& f : faces_) {
    if (!f->IsInBoundingBox(global)) {
      continue;
    }
    auto [dist, param] = f->Project(global);
    if (dist <= 1e-5 && f->IsInside(param)) {
      result.emplace_back(f, std::move(param));
    }
  }
  return result;
}

std::vector<
    std::pair<std::shared_ptr<const interface::BrepSurface>, Eigen::Matrix2Xd>>
OcctBrepModel::FindSurfacesMulti(const Eigen::Matrix3Xd& global) const {
  std::vector<std::pair<std::shared_ptr<const interface::BrepSurface>,
                        Eigen::Matrix2Xd>>
      result;
  Eigen::Matrix2Xd result_coord;
  for (const auto& f : faces_) {
    // 1) Check if all points are inside the bounding box:
    if ([&]() {
          for (int i = 0; i < global.cols(); ++i) {
            if (!f->IsInBoundingBox(global.col(i))) {
              return true;
            }
          }
          return false;
        }()) {
      continue;
    }

    // project points onto face and check if it is within bounds:
    for (unsigned i = 0; i < global.cols(); ++i) {
      auto [dist, param] = f->Project(global.col(i));
      if (dist > 1e-5 || !f->IsInside(param)) {
        break;
      }
      if (result_coord.size() == 0) {
        result_coord.resize(2, global.cols());
      }
      result_coord.col(i) = param;
      if (i + 1 == global.cols()) {
        result.emplace_back(f, std::move(result_coord));
        result_coord = Eigen::Matrix2Xd();
      }
    }
  }
  return result;
}

}  // namespace lf::brep::occt

/**
 * @file
 * @brief Defines the OcctBrepModel class
 * @author Raffael Casagrande
 * @date   2020-11-06 05:18:16
 * @copyright MIT License
 */

#ifndef __c22f41716f574518a166a992467c0172
#define __c22f41716f574518a166a992467c0172

#include <Bnd_OBB.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include "lf/brep/interface/interface.h"
#include "occt_brep_curve.h"
#include "occt_brep_surface.h"

namespace lf::brep::occt {

class OcctBrepModel : public interface::BrepModel {
 public:
  explicit OcctBrepModel(std::string_view filename);

  /**
   * @brief Find all Curves that go through the given point.
   */
  [[nodiscard]] std::vector<
      std::pair<std::shared_ptr<const interface::BrepCurve>, double>>
  FindCurves(const Eigen::Vector3d &global) const override;

  /**
   * @brief Find all Curves that go through *all* the given points (as column
   * vectors).
   */
  [[nodiscard]] std::vector<std::pair<
      std::shared_ptr<const interface::BrepCurve>, Eigen::RowVectorXd>>
  FindCurvesMulti(const Eigen::Matrix3Xd &global) const override;

  /**
   * @brief Find all Surfaces that contain the given point.
   */
  [[nodiscard]] std::vector<
      std::pair<std::shared_ptr<const interface::BrepSurface>, Eigen::Vector2d>>
  FindSurfaces(const Eigen::Vector3d &global) const override;

  /**
   * @brief Find all Surfaces that contain *all* the given points (as column
   * vectors)
   */
  [[nodiscard]] std::vector<std::pair<
      std::shared_ptr<const interface::BrepSurface>, Eigen::Matrix2Xd>>
  FindSurfacesMulti(const Eigen::Matrix3Xd &global) const override;

  /**
   * @brief Total number of Curves in the model
   */
  [[nodiscard]] base::size_type NumCurves() const override {
    return edges_.size();
  }

  /**
   * @brief Total number of Surfaces in the model.
   */
  [[nodiscard]] base::size_type NumSurfaces() const override {
    return faces_.size();
  }

 private:
  TopoDS_Shape shape_;
  std::vector<std::shared_ptr<OcctBrepCurve>> edges_;
  std::vector<std::shared_ptr<OcctBrepSurface>> faces_;
};

}  // namespace lf::brep::occt

#endif  // __c22f41716f574518a166a992467c0172

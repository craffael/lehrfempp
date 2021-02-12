/**
 * @file
 * @brief Defines the OCCTSurfaceGeometry class
 * @author Raffael Casagrande
 * @date   2020-11-20 03:15:45
 * @copyright MIT License
 */

#ifndef __9b2aa16d069a447eaa411ed6c262f4f1
#define __9b2aa16d069a447eaa411ed6c262f4f1

#include <lf/brep/interface/interface.h>

#include <Bnd_OBB.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Face.hxx>

namespace lf::brep::occt {

class OcctBrepSurface : public interface::BrepGeometry {
 public:
  OcctBrepSurface(TopoDS_Face&& face);
  OcctBrepSurface(const OcctBrepSurface&) = default;
  OcctBrepSurface(OcctBrepSurface&&) = default;
  OcctBrepSurface& operator=(const OcctBrepSurface&) = default;
  OcctBrepSurface& operator=(OcctBrepSurface&&) = default;
  ~OcctBrepSurface() = default;

  [[nodiscard]] base::dim_t DimGlobal() const override { return 3; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 2; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] std::pair<double, Eigen::VectorXd> Project(
      const Eigen::VectorXd& global) const override;
  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const override;

  [[nodiscard]] bool IsInside(const Eigen::VectorXd& local) const override;

  // OCCT specific member functions:

  [[nodiscard]] const TopoDS_Face& Face() const { return face_; }

  [[nodiscard]] Eigen::VectorXd Periods() const override {
    Eigen::VectorXd result(2);
    if (surface_->IsUPeriodic()) {
      result(0) = surface_->UPeriod();
    } else {
      result(0) = 0;
    }
    if (surface_->IsVPeriodic()) {
      result(0) = surface_->VPeriod();
    } else {
      result(0) = 0;
    }
    return result;
  }

 private:
  TopoDS_Face face_;
  Bnd_OBB obb_;
  opencascade::handle<Geom_Surface> surface_;
};

}  // namespace lf::brep::occt

#endif  // __9b2aa16d069a447eaa411ed6c262f4f1

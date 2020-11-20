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

class OcctSurfaceGeometry : public interface::BrepGeometry {
 public:
  OcctSurfaceGeometry(TopoDS_Face&& face);
  OcctSurfaceGeometry(const OcctSurfaceGeometry&) = default;
  OcctSurfaceGeometry(OcctSurfaceGeometry&&) = default;
  OcctSurfaceGeometry& operator=(const OcctSurfaceGeometry&) = default;
  OcctSurfaceGeometry& operator=(OcctSurfaceGeometry&&) = default;
  ~OcctSurfaceGeometry() = default;

  [[nodiscard]] base::dim_t DimGlobal() const override { return 3; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 2; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] std::pair<Eigen::VectorXd, Eigen::MatrixXd> Project(
      const Eigen::MatrixXd& global) const override;
  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const override;
  [[nodiscard]] std::vector<bool> IsInside(
      const Eigen::MatrixXd& local) const override;

 private:
  TopoDS_Face face_;
  Bnd_OBB obb_;
  opencascade::handle<Geom_Surface> surface_;
};

}  // namespace lf::brep::occt

#endif  // __9b2aa16d069a447eaa411ed6c262f4f1

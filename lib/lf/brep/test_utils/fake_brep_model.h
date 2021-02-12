/**
 * @file
 * @brief Implementation of a interface::BrepModel for testing purposes.
 * @author Raffael Casagrande
 * @date   2021-02-11 02:06:25
 * @copyright MIT License
 */

#ifndef __de1f797367a54114966ff8a962d027c4
#define __de1f797367a54114966ff8a962d027c4

#include <lf/brep/interface/interface.h>

namespace lf::brep::test_utils {

class FakeBrepModel : public interface::BrepModel {
 public:
  [[nodiscard]] std::vector<
      std::pair<std::shared_ptr<const interface::BrepGeometry>, double>>
  FindCurves(const Eigen::Vector3d& global) const override;

  [[nodiscard]] std::vector<std::pair<
      std::shared_ptr<const interface::BrepGeometry>, Eigen::RowVectorXd>>
  FindCurvesMulti(const Eigen::Matrix3Xd& global) const override;

  [[nodiscard]] std::vector<std::pair<
      std::shared_ptr<const interface::BrepGeometry>, Eigen::Vector2d>>
  FindSurfaces(const Eigen::Vector3d& global) const override;

  [[nodiscard]] std::vector<std::pair<
      std::shared_ptr<const interface::BrepGeometry>, Eigen::Matrix2Xd>>
  FindSurfacesMulti(const Eigen::Matrix3Xd& global) const override;

  [[nodiscard]] base::size_type NumCurves() const override {
    return curves_.size();
  }
  [[nodiscard]] base::size_type NumSurfaces() const override {
    return surfaces_.size();
  }

  void AddCurve(std::shared_ptr<const interface::BrepGeometry> curve) {
    LF_ASSERT_MSG(curve->DimLocal() == 1, "This is not a curve.");
    curves_.push_back(std::move(curve));
  }

  void AddSurface(std::shared_ptr<const interface::BrepGeometry> surface) {
    LF_ASSERT_MSG(surface->DimLocal() == 2, "This is not a surface.");
    surfaces_.push_back(std::move(surface));
  }

 private:
  std::vector<std::shared_ptr<const interface::BrepGeometry>> curves_;
  std::vector<std::shared_ptr<const interface::BrepGeometry>> surfaces_;

  template <class COORD>
  auto static FindGeometry(
      const std::vector<std::shared_ptr<const interface::BrepGeometry>>&
          geometries,
      const Eigen::Vector3d& global);

  template <class COORD>
  auto FindGeometryMulti(
      const std::vector<std::shared_ptr<const interface::BrepGeometry>>&
          geometries,
      const Eigen::Matrix3Xd& global) const;
};

}  // namespace lf::brep::test_utils

#endif  // __de1f797367a54114966ff8a962d027c4

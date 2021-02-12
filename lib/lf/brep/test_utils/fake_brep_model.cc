/**
 * @file
 * @brief Implementation of FakeBrepModel.
 * @author Raffael Casagrande
 * @date   2021-02-11 02:09:30
 * @copyright MIT License
 */

#include "fake_brep_model.h"

namespace lf::brep::test_utils {

template <class COORD>
auto FakeBrepModel::FindGeometry(
    const std::vector<std::shared_ptr<const interface::BrepGeometry>>&
        geometries,
    const Eigen::Vector3d& global) {
  std::vector<std::pair<std::shared_ptr<const interface::BrepGeometry>, COORD>>
      result;
  for (const auto& c : geometries) {
    auto [dist, param] = c->Project(global);
    if (dist < 1e-6 && c->IsInside(param)) {
      if constexpr (std::is_same_v<COORD, double>) {
        result.emplace_back(c, param(0));
      } else {
        result.emplace_back(c, param);
      }
    }
  }
  return result;
}

template <class COORD>
auto FakeBrepModel::FindGeometryMulti(
    const std::vector<std::shared_ptr<const interface::BrepGeometry>>&
        geometries,
    const Eigen::Matrix3Xd& global) const {
  std::vector<std::pair<std::shared_ptr<const interface::BrepGeometry>, COORD>>
      result;
  Eigen::MatrixXd param(geometries == curves_ ? 1 : 2, global.cols());
  for (const auto& c : geometries) {
    bool found_all = true;
    for (int i = 0; i < global.cols(); ++i) {
      auto [dist, p] = c->Project(global.col(i));
      if (dist > 1e-6 || c->IsInside(p)) {
        found_all = false;
        break;
      }
      param.col(i) = p;
    }
    if (found_all) {
      result.emplace_back(c, param);
    }
  }
  return result;
}

std::vector<std::pair<std::shared_ptr<const interface::BrepGeometry>, double>>
FakeBrepModel::FindCurves(const Eigen::Vector3d& global) const {
  return FindGeometry<double>(curves_, global);
}

std::vector<std::pair<std::shared_ptr<const interface::BrepGeometry>,
                      Eigen::RowVectorXd>>
FakeBrepModel::FindCurvesMulti(const Eigen::Matrix3Xd& global) const {
  return FindGeometryMulti<Eigen::RowVectorXd>(curves_, global);
}

std::vector<
    std::pair<std::shared_ptr<const interface::BrepGeometry>, Eigen::Vector2d>>
FakeBrepModel::FindSurfaces(const Eigen::Vector3d& global) const {
  return FindGeometry<Eigen::Vector2d>(surfaces_, global);
}

std::vector<
    std::pair<std::shared_ptr<const interface::BrepGeometry>, Eigen::Matrix2Xd>>
FakeBrepModel::FindSurfacesMulti(const Eigen::Matrix3Xd& global) const {
  return FindGeometryMulti<Eigen::Matrix2Xd>(surfaces_, global);
}
}  // namespace lf::brep::test_utils

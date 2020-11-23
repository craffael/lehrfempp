/**
 * @file
 * @brief Some utility functions to test the lf::brep::occt module
 * @author Raffael Casagrande
 * @date   2020-11-11 02:50:52
 * @copyright MIT License
 */

#ifndef __15c6adac3fc54dd8b4b887262e606ea6
#define __15c6adac3fc54dd8b4b887262e606ea6

#include "lf/brep/occt/occt.h"

#include <filesystem>

namespace lf::brep::occt::test {

inline std::unique_ptr<OcctBrepModel> LoadModel(std::string_view filename) {
  std::filesystem::path here = __FILE__;
  return std::make_unique<OcctBrepModel>(
      (here.remove_filename() / "brep_models" / filename).string());
}

template<class F>
Eigen::Matrix3Xd approxJacobian(const F& f, Eigen::VectorXd local) {
  Eigen::Matrix3Xd result(3, local.rows());
  static const double eps = 1e-5; // See Timothy Sauter, Numerical Analysis
  for (int i = 0; i < local.rows(); ++i) {
    Eigen::VectorXd dh=
        eps * Eigen::MatrixXd::Identity(local.rows(), local.rows()).col(i);
    result.col(i) = (f(local + dh) - f(local - dh)) / (2 * eps);
  }
  return result;
}

inline void CheckGeometry(interface::BrepGeometry const* geom, Eigen::MatrixXd local) {
  EXPECT_EQ(geom->DimGlobal(), 3);
  auto global = geom->GlobalMulti(local);
  EXPECT_EQ(global.cols(), local.cols());
  EXPECT_EQ(global.rows(), 3);

  auto bb = geom->IsInBoundingBoxMulti(global);
  for (auto c : bb) {
    EXPECT_TRUE(c);
  }

  auto jacobian = geom->JacobianMulti(local);

  if (auto curve = dynamic_cast<interface::BrepCurve const*>(geom); curve) {
    EXPECT_EQ(curve->DimLocal(), 1);
    for (int i = 0; i < local.cols(); ++i) {
      EXPECT_TRUE(curve->IsInBoundingBoxSingle(global.col(i)));
      EXPECT_TRUE(curve->IsInside(local(0,i)));
      EXPECT_TRUE(curve->GlobalSingle(local(0,i)).isApprox(global.col(i)));
      auto jacApprox = approxJacobian(
          [&](const auto& x) { return curve->GlobalSingle(x(0, 0)); },
          local.col(i));
      ASSERT_TRUE(
          curve->JacobianSingle(local(0,i)).isApprox(jacobian.col(i)));
      ASSERT_TRUE(jacApprox.isApprox(jacobian.col(i), 1e-5));

      auto [dist, pl] = curve->Project(global.col(i));
      ASSERT_TRUE(std::abs(dist) < 1e-6);
      ASSERT_TRUE(std::abs(pl - local(0, i)) < 1e-6);
    }
  } else {
    auto surf = dynamic_cast<interface::BrepSurface const*>(geom);
    EXPECT_EQ(surf->DimLocal(), 2);
    EXPECT_TRUE(surf);
    for (int i=0; i<local.cols(); ++i) {
      EXPECT_TRUE(surf->IsInBoundingBox(global.col(i)));
      EXPECT_TRUE(surf->IsInside(local.col(i)));
      EXPECT_TRUE(surf->Global(local.col(i)).isApprox(global.col(i)));

      auto jacApprox = approxJacobian(
          [&](const auto& x) { return surf->Global(x); }, local.col(i));
      EXPECT_TRUE(jacApprox.isApprox(jacobian.block<3,2>(0,2*i), 1e-5));
      EXPECT_TRUE(surf->Jacobian(local.col(i))
                      .isApprox(jacobian.block<3, 2>(0, 2 * i)));

      auto [dist, pl] = surf->Project(global.col(i));
      ASSERT_TRUE(std::abs(dist) < 1e-6);
      ASSERT_TRUE(pl.isApprox(local.col(i)));
    }
  }
}

}  // namespace lf::brep::occt::test

#endif  // __15c6adac3fc54dd8b4b887262e606ea6

/**
 * @file
 * @brief implementation from check_brep_geometry.h
 * @author Raffael Casagrande
 * @date   2021-01-29 10:22:23
 * @copyright MIT License
 */

#include "check_brep_geometry.h"

namespace lf::brep::test_utils {
void CheckBrepGeometry(interface::BrepGeometry const& geom,
                       Eigen::MatrixXd local) {
  EXPECT_EQ(geom.DimGlobal(), 3);
  auto global = geom.GlobalMulti(local);
  EXPECT_EQ(global.cols(), local.cols());
  EXPECT_EQ(global.rows(), 3);

  auto bb = geom.IsInBoundingBoxMulti(global);
  for (auto c : bb) {
    EXPECT_TRUE(c);
  }

  auto jacobian = geom.JacobianMulti(local);

  if (auto curve = dynamic_cast<interface::BrepCurve const*>(&geom); curve) {
    EXPECT_EQ(curve->DimLocal(), 1);
    for (int i = 0; i < local.cols(); ++i) {
      EXPECT_TRUE(curve->IsInBoundingBoxSingle(global.col(i)));
      EXPECT_TRUE(curve->IsInside(local(0, i)));
      EXPECT_TRUE(curve->GlobalSingle(local(0, i)).isApprox(global.col(i)));
      auto jacApprox = approxJacobian(
          [&](const auto& x) { return curve->GlobalSingle(x(0, 0)); },
          local.col(i));
      ASSERT_TRUE(curve->JacobianSingle(local(0, i)).isApprox(jacobian.col(i)));
      ASSERT_TRUE(jacApprox.isApprox(jacobian.col(i), 1e-5));

      auto [dist, pl] = curve->Project(global.col(i));
      ASSERT_TRUE(std::abs(dist) < 1e-6);
      ASSERT_TRUE(std::abs(pl - local(0, i)) < 1e-6);
    }
  } else {
    auto surf = dynamic_cast<interface::BrepSurface const*>(&geom);
    EXPECT_EQ(surf->DimLocal(), 2);
    EXPECT_TRUE(surf);
    for (int i = 0; i < local.cols(); ++i) {
      EXPECT_TRUE(surf->IsInBoundingBox(global.col(i)));
      EXPECT_TRUE(surf->IsInside(local.col(i)));
      EXPECT_TRUE(surf->GlobalSingle(local.col(i)).isApprox(global.col(i)));

      auto jacApprox = approxJacobian(
          [&](const auto& x) { return surf->GlobalSingle(x); }, local.col(i));
      EXPECT_TRUE(jacApprox.isApprox(jacobian.block<3, 2>(0, 2 * i), 1e-5));
      EXPECT_TRUE(surf->JacobianSingle(local.col(i))
                      .isApprox(jacobian.block<3, 2>(0, 2 * i)));

      auto [dist, pl] = surf->Project(global.col(i));
      ASSERT_TRUE(std::abs(dist) < 1e-6);
      ASSERT_TRUE(pl.isApprox(local.col(i)));
    }
  }
}
}  // namespace lf::brep::test_utils

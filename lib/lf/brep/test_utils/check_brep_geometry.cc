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
  auto global = geom.Global(local);
  EXPECT_EQ(global.cols(), local.cols());
  EXPECT_EQ(global.rows(), 3);

  auto bb = geom.IsInBoundingBox(global);
  for (auto c : bb) {
    EXPECT_TRUE(c);
  }

  auto jacobian = geom.Jacobian(local);

  if (geom.DimLocal() == 1) {
    for (int i = 0; i < local.cols(); ++i) {
      EXPECT_TRUE(geom.IsInside(local.col(i)));

      auto jacApprox = approxJacobian(
          [&](const auto& x) { return geom.Global(x); }, local.col(i));
      ASSERT_TRUE(jacApprox.isApprox(jacobian.col(i), 1e-5));

      auto [dist, pl] = geom.Project(global.col(i));
      ASSERT_EQ(pl.rows(), 1);
      ASSERT_TRUE(std::abs(dist) < 1e-6);
      ASSERT_TRUE(std::abs(pl(0, 0) - local(0, i)) < 1e-6);
    }
  } else {
    EXPECT_EQ(geom.DimLocal(), 2);
    for (int i = 0; i < local.cols(); ++i) {
      EXPECT_TRUE(geom.IsInside(local.col(i)));

      auto jacApprox = approxJacobian(
          [&](const auto& x) { return geom.Global(x); }, local.col(i));
      EXPECT_TRUE(jacApprox.isApprox(jacobian.block<3, 2>(0, 2 * i), 1e-5));

      auto [dist, pl] = geom.Project(global.col(i));
      ASSERT_EQ(pl.rows(), 2);
      ASSERT_TRUE(std::abs(dist) < 1e-6);
      ASSERT_TRUE(pl.isApprox(local.col(i)));
    }
  }

  // Check periodicity:
  auto periods = geom.Periods();
  ASSERT_EQ(periods.rows(), geom.DimLocal());
  for (int d = 0; d < geom.DimLocal(); ++d) {
    if (periods[d] != 0) {
      auto unit =
          (Eigen::VectorXd::Unit(periods.rows(), d) * periods[d]).eval();
      ASSERT_TRUE(global.isApprox(
          geom.Global(local + unit.replicate(1, local.cols()))));
    }
  }
}
}  // namespace lf::brep::test_utils

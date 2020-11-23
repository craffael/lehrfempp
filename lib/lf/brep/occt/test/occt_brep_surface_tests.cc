/** 
 * @file 
 * @brief Test implementation of OcctBrepSurface
 * @author Raffael Casagrande
 * @date   2020-11-23 01:04:03
 * @copyright MIT License
 */

#include <gtest/gtest.h>

#include <lf/brep/occt/occt.h>

#include "utils.h"

namespace lf::brep::occt::test {

TEST(occt, cubeSurfaceProjection) {
  auto model = LoadModel("cube.brep");

  Eigen::MatrixXd global(3, 3);
  global.col(0) = Eigen::Vector3d(10, 10, 10);
  global.col(1) = Eigen::Vector3d(10, 10, 0);
  global.col(2) = Eigen::Vector3d(0, 10, 0);
  auto surf_result = model->FindSurfacesMulti(global);
  ASSERT_EQ(surf_result.size(), 1);
  auto surface = surf_result[0].first;

  // project the point (5,10,5):
  Eigen::Vector3d p(5, 10, 5);
  auto [dist, param] = surface->Project(p);
  EXPECT_LT(dist, 1e-6);
  EXPECT_TRUE(surface->Global(param).isApprox(p));

  // Project the point (5,15,5):
  p = Eigen::Vector3d(5, 15, 5);
  std::tie(dist, param) = surface->Project(p);
  EXPECT_LT(std::abs(dist - 5), 1e-5);
  EXPECT_TRUE(surface->Global(param).isApprox(Eigen::Vector3d(5,10,5)));

  // project a point that lies slightly outside:
  p = Eigen::Vector3d(11, 10, 0);
  std::tie(dist, param) = surface->Project(p);
  EXPECT_LT(dist, 1e-5);
  EXPECT_FALSE(surface->IsInside(param));
}

TEST(OCCT, hollowHemiespherSurfProjection) {
  auto model = LoadModel("hollow_hemisphere.brep");

  // find the outer shell:
  auto surface_result = model->FindSurfaces(Eigen::Vector3d(0, 0, -5));
  ASSERT_EQ(surface_result.size(), 1);
  auto surface = surface_result[0].first;

  // Project point (-6,-6,-6) onto the surface:
  Eigen::Vector3d p(-6, -6, -6);
  auto [dist, param] = surface->Project(p);
  EXPECT_GT(dist, 1);
  EXPECT_TRUE(surface->Global(param).isApprox(p / p.norm() * 5));

  // Project point (6,6,6) onto the surface:
  // -> Here it's not clear if a point within the parameter bounds will be returned or not.
  p = Eigen::Vector3d(6, 6, 6);
  std::tie(dist, param) = surface->Project(p);
  //std::cout << surface->Global(param) << std::endl;
  //EXPECT_TRUE(surface->IsInside(param));


}
}





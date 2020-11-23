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

inline void CheckGeometry(interface::BrepGeometry const* geom, Eigen::MatrixXd local) {
  auto global = geom->Global(local);
  EXPECT_EQ(global.cols(), local.cols());
  EXPECT_EQ(global.rows(), 3);

  auto bb = geom->IsInBoundingBox(global);
  for (auto c : bb) {
    EXPECT_TRUE(c);
  }

  if (auto curve = dynamic_cast<interface::BrepCurve const*>(geom); curve) {
    for (int i = 0; i < local.cols(); ++i) {
      EXPECT_TRUE(curve->IsInBoundingBoxSingle(global.col(i)));
      EXPECT_TRUE(curve->IsInside(local(0,i)));
      EXPECT_TRUE(curve->GlobalSingle(local(0,i)).isApprox(global.col(i)));
    }
  } else {
    auto surf = dynamic_cast<interface::BrepSurface const*>(geom);
    EXPECT_TRUE(surf);
    for (int i=0; i<local.cols(); ++i) {
      EXPECT_TRUE(surf->IsInBoundingBoxSingle(global.col(i)));
      EXPECT_TRUE(surf->IsInside(local.col(i)));
      EXPECT_TRUE(surf->GlobalSingle(local.col(i)).isApprox(global.col(i)));
    }
  }
  
  
}

}  // namespace lf::brep::occt::test

#endif  // __15c6adac3fc54dd8b4b887262e606ea6

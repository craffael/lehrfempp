#include <gtest/gtest.h>

#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_function_transfer.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <cmath>
#include <functional>
#include <vector>

TEST(lf_refinement, MeshFunctionTransferScalar) {
  // Declare some arbitrary functions from R^2 to R that will be transferred
  std::vector<std::function<double(const Eigen::Vector2d&)>> functions;
  functions.push_back([](const Eigen::Vector2d& x) { return 1; });
  functions.push_back([](const Eigen::Vector2d& x) { return x[0]; });
  functions.push_back([](const Eigen::Vector2d& x) { return x[1]; });
  functions.push_back([](const Eigen::Vector2d& x) { return x[0] + x[1]; });
  functions.push_back(
      [](const Eigen::Vector2d& x) { return std::sin(x[0] + x[1]); });
  // Generate a test mesh and put that into a refinement hierarchy
  for (int selector = 0; selector < 9; ++selector) {
    // Construct a MeshHierarchy from the test mesh
    const auto mesh_coarse =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    auto hierarchy =
        lf::refinement::MeshHierarchy(mesh_coarse, std::move(factory));
    hierarchy.RefineRegular();
    // Transfer all test functions
    for (const auto& f : functions) {
      lf::mesh::utils::MeshFunctionGlobal mf_global(f);
      // Transfer the test function
      lf::refinement::MeshFunctionTransfer mf_fine(hierarchy, mf_global, 0);
      // Compare the two mesh functions at the midpoints of the fine mesh
      for (const auto cell : hierarchy.getMesh(1)->Entities(0)) {
        const Eigen::Vector2d midpoint = Eigen::Vector2d::Constant(1. / 3);
        const double value_global = mf_global(*cell, midpoint)[0];
        const double value_transferred = mf_fine(*cell, midpoint)[0];
        ASSERT_NEAR(value_global, value_transferred, 1e-10)
            << "value_global=" << value_global
            << "\tvalue_transferred=" << value_transferred;
      }
    }
  }
}

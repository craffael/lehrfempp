/**
 * @file
 * @brief Check whether MeshFunctionFE works the way it should.
 * @author Raffael Casagrande
 * @date   2019-01-19 05:58:08
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/test/mesh_function_utils.h"

namespace lf::uscalfe::test {

TEST(lf_fe_meshFunctionFE, Projection) {
  // This test projects a linear mesh function onto a fe space and compares
  // the MeshFunctionFE with the original mesh function.

  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  auto mf_linear = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] + 2 * x[1]; });

  auto fespaceO1 = std::make_shared<FeSpaceLagrangeO1<double>>(mesh);

  auto projected = lf::fe::NodalProjection(*fespaceO1, mf_linear);
  auto mf_projected = lf::fe::MeshFunctionFE(fespaceO1, projected);
  mesh::utils::test::checkMeshFunctionEqual(*mesh, mf_linear, mf_projected);
}

TEST(lf_fe_meshFunctionFE, ProjectionQuad) {
  // This test projects a quadratic mesh function onto a fe space and compares
  // the MeshFunctionFE with the original mesh function.

  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // A quadratic polynomial that can exactly be represented in the FE space
  auto mf_quad = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] * x[1] + 2 * x[1]; });

  auto fespaceO2 = std::make_shared<FeSpaceLagrangeO2<double>>(mesh);

  auto projected = NodalProjection(*fespaceO2, mf_quad);
  auto mf_projected = lf::fe::MeshFunctionFE(fespaceO2, projected);
  mesh::utils::test::checkMeshFunctionEqual(*mesh, mf_quad, mf_projected);
}

TEST(lf_fe_meshFunctionFE, ProjectionCubic) {
  // This test projects a quadratic mesh function onto a fe space and compares
  // the MeshFunctionFE with the original mesh function.

  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // A cubic polynomial that can exactly be represented in the FE space
  auto mf_cubic = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
    return x[0] * x[1] * x[1] + x[0] * x[0] + 2 * x[1];
  });

  auto fespaceO3 = std::make_shared<FeSpaceLagrangeO3<double>>(mesh);

  auto projected = NodalProjection(*fespaceO3, mf_cubic);
  auto mf_projected = lf::fe::MeshFunctionFE(fespaceO3, projected);
  mesh::utils::test::checkMeshFunctionEqual(*mesh, mf_cubic, mf_projected);
}

}  // namespace lf::uscalfe::test

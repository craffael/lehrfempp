/**
 * @file
 * @brief Check MeshFunctionGradFE
 * @author Raffael Casagrande
 * @date   2019-01-19 06:50:42
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/uscalfe/uscalfe.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/test/mesh_function_utils.h"

namespace lf::uscalfe::test {

TEST(lf_fe_meshFunctionGradFE, ProjectionTest) {
  // This test projects a linear mesh function onto a fe space and compares
  // the MeshFunctionFE with the original mesh function.

  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  auto mf_linear = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] + 2 * x[1]; });
  auto mf_grad =
      mesh::utils::MeshFunctionConstant(Eigen::VectorXd(Eigen::Vector2d(1, 2)));

  auto fespaceO1 = std::make_shared<FeSpaceLagrangeO1<double>>(mesh);

  auto projected = lf::fe::NodalProjection(*fespaceO1, mf_linear);
  auto mf_grad_projected = lf::fe::MeshFunctionGradFE(fespaceO1, projected);
  mesh::utils::test::checkMeshFunctionEqual(*mesh, mf_grad, mf_grad_projected);
}

TEST(lf_fe_meshFunctionGradFE, ProjectionTestQuad) {
  // This test projects a linear mesh function onto a quadratic Lagrangian FE
  // space and compares the MeshFunctionFE with the original mesh function.

  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // A two-variate polynomials of degree 2
  // This function can exactly be represented in the FE space
  auto mf_quad = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] * x[1] + 2 * x[1]; });
  auto mf_grad = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
    return Eigen::VectorXd(Eigen::Vector2d(x[1], x[0] + 2.0));
  });

  auto fespaceO2 = std::make_shared<FeSpaceLagrangeO2<double>>(mesh);

  auto projected = NodalProjection(*fespaceO2, mf_quad);
  auto mf_grad_projected = lf::fe::MeshFunctionGradFE(fespaceO2, projected);
  mesh::utils::test::checkMeshFunctionEqual(*mesh, mf_grad, mf_grad_projected);
}

}  // namespace lf::uscalfe::test

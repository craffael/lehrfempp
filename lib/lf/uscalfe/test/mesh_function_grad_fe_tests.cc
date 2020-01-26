/**
 * @file
 * @brief Check MeshFunctionGradFE
 * @author Raffael Casagrande
 * @date   2019-01-19 06:50:42
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/uscalfe/uscalfe.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/test/mesh_function_utils.h"

namespace lf::uscalfe::test {

TEST(meshFunctionGradFE, ProjectionTest) {
  // This test projects a linear mesh function onto a fe space and compares
  // the MeshFunctionFE with the original mesh function.

  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  auto mf_linear = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] + 2 * x[1]; });
  auto mf_grad =
      mesh::utils::MeshFunctionConstant(Eigen::VectorXd(Eigen::Vector2d(1, 2)));

  auto fespaceO1 = std::make_shared<FeSpaceLagrangeO1<double>>(mesh);

  auto projected = NodalProjection(*fespaceO1, mf_linear);
  auto mf_grad_projected = MeshFunctionGradFE(fespaceO1, projected);
  mesh::utils::test::checkMeshFunctionEqual(*mesh, mf_grad, mf_grad_projected);
}

}  // namespace lf::uscalfe::test

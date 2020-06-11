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

TEST(meshFunctionFE, Projection) {
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

}  // namespace lf::uscalfe::test

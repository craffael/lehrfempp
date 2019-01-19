/**
 * @file
 * @brief Check whether MeshFunctionFE works the way it should.
 * @author Raffael Casagrande
 * @date   2019-01-19 05:58:08
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/uscalfe/uscalfe.h>
#include "mesh_function_utils.h"

namespace lf::uscalfe::test {

TEST(meshFunctionFE, Projection) {
  // This test projects a linear mesh function onto a fe space and compares
  // the MeshFunctionFE with the original mesh function.

  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();
  auto mf_linear = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x[0] + 2 * x[1]; });

  auto fespaceO1 = std::make_shared<FeSpaceLagrangeO1<double>>(mesh);

  auto projected = NodalProjection(*fespaceO1, mf_linear);
  auto mf_projected = MeshFunctionFE(fespaceO1, projected);
  checkMeshFunctionEqual(*mesh, mf_linear, mf_projected);
}

}  // namespace lf::uscalfe::test

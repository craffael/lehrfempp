/**
 * @file
 * @brief Check that the classes in fe_tools.h work as expected.
 * @author Raffael Casagrande
 * @date   2019-01-19 09:18:35
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::uscalfe::test {

TEST(feTools, IntegrateMeshFunction) {
  // [0,3]^2
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  io::VtkWriter vtk_writer(mesh, "mesh.vtk");

  // scalar valued mesh function
  auto mfScalar = MeshFunctionGlobal([](const Eigen::Vector2d& x) {
    return std::sin(base::kPi * x[0]) * x[1];
  });
  auto intScalar = IntegrateMeshFunction(*mesh, mfScalar, 20);
  EXPECT_FLOAT_EQ(intScalar, 9 / base::kPi);

  // vector valued mesh function
  auto mfVec = MeshFunctionGlobal([](const Eigen::Vector2d& x) {
    return Eigen::Vector3d(x[0], x[0] * x[1], std::cos(base::kPi * x[0]));
  });
  auto intVec = IntegrateMeshFunction(*mesh, mfVec, 20);
  EXPECT_FLOAT_EQ(intVec[0], 27. / 2.);
  EXPECT_FLOAT_EQ(intVec[1], 81. / 4.);
  EXPECT_LT(intVec[2], 1e-10);

  // dynamic matrix sized mesh function
  auto mfMatrixDyn =
      MeshFunctionGlobal([](const Eigen::Vector2d& x) -> Eigen::MatrixXd {
        return (Eigen::MatrixXd(2, 3) << 1, x[0], x[1], x[0] * x[1],
                x[0] * x[0], x[1] * x[1])
            .finished();
      });
  auto intMatrixDyn = IntegrateMeshFunction(*mesh, mfMatrixDyn, 10);
  EXPECT_FLOAT_EQ(intMatrixDyn(0, 0), 9.);
  EXPECT_FLOAT_EQ(intMatrixDyn(0, 1), 27. / 2.);
  EXPECT_FLOAT_EQ(intMatrixDyn(0, 2), 27. / 2.);
  EXPECT_FLOAT_EQ(intMatrixDyn(1, 0), 81. / 4.);
  EXPECT_FLOAT_EQ(intMatrixDyn(1, 1), 27.);
  EXPECT_FLOAT_EQ(intMatrixDyn(1, 2), 27.);
}

}  // namespace lf::uscalfe::test

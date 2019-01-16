/**
 * @file
 * @brief Check if all binary mesh functions work as expected.
 * @author Raffael Casagrande
 * @date   2019-01-13 05:28:53
 * @copyright MIT License
 */

#include <gtest/gtest.h>

#include <lf/io/test_utils/read_mesh.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>
#include "mesh_function_utils.h"

namespace lf::uscalfe::test {

class X {
 public:
  X() {}
  X(int x) : x_(x) {}

  X operator+(const X& rhs) const { return X(x_ + rhs.x_); }
  X operator-(const X& rhs) const { return X(x_ - rhs.x_); }
  bool operator==(const X& rhs) const { return x_ == rhs.x_; }

 private:
  int x_;
};

auto mfA = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) { return x[0] * x[0] + x[1]; });
auto mfB =
    MeshFunctionGlobal([](const Eigen::Vector2d& x) { return x[0] + x[1]; });
auto mfVectorA = MeshFunctionGlobal([](const Eigen::Vector2d& x) { return x; });
auto mfVectorB = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) { return Eigen::Vector2d(x[0], 2 * x[0]); });

auto mfMatrixA =
    MeshFunctionGlobal([](auto x) { return (x * x.transpose()).eval(); });
auto mfMatrixB = MeshFunctionGlobal([](auto x) -> Eigen::Matrix2d {
  return (Eigen::Matrix2d() << 0, 1, x[0], x[1]).finished();
});

auto mfXA = MeshFunctionConstant(X(1));
auto mfXB = MeshFunctionConstant(X(2));

TEST(meshFunctionBinary, Addition) {
  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();

  auto sum = mfA + mfB;
  checkMeshFunctionEqual(*mesh, sum, MeshFunctionGlobal([](auto x) {
    return x[0] * x[0] + 2 * x[1] + x[0];
  }));

  auto vectorSum = mfVectorA + mfVectorB;
  checkMeshFunctionEqual(*mesh, vectorSum, MeshFunctionGlobal([](auto x) {
    return Eigen::Vector2d(2 * x[0], x[1] + 2 * x[0]);
  }));

  auto matrixSum = mfMatrixA + mfMatrixB;
  checkMeshFunctionEqual(*mesh, matrixSum, MeshFunctionGlobal([](auto x) {
    return (Eigen::Matrix2d() << x[0] * x[0], 1 + x[0] * x[1],
            x[0] * x[1] + x[0], x[1] * x[1] + x[1])
        .finished();
  }));

  auto xSum = mfXA + mfXB;
  checkMeshFunctionEqual(*mesh, xSum, MeshFunctionConstant(X(3)));
}

TEST(meshFunctionBinary, Subtraction) {
  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();

  auto sub = mfA - mfB;
  checkMeshFunctionEqual(*mesh, sub, MeshFunctionGlobal([](auto x) {
    return x[0] * x[0] - x[0];
  }));

  auto subVector = mfVectorA - mfVectorB;
  checkMeshFunctionEqual(*mesh, subVector, MeshFunctionGlobal([](auto x) {
    return Eigen::Vector2d(0., x[1] - 2 * x[0]);
  }));

  auto subMatrix = mfMatrixA - mfMatrixB;
  checkMeshFunctionEqual(*mesh, subMatrix, MeshFunctionGlobal([](auto x) {
    return (Eigen::Matrix2d() << x[0] * x[0], x[0] * x[1] - 1,
            x[0] * x[1] - x[0], x[1] * x[1] - x[1])
        .finished();
  }));

  auto xSub = mfXA - mfXB;
  checkMeshFunctionEqual(*mesh, xSub, MeshFunctionConstant(X(-1)));
}

}  // namespace lf::uscalfe::test

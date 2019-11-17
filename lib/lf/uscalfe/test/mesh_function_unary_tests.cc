/**
 * @file
 * @brief Test the unary mesh functions
 * @author Raffael Casagrande
 * @date   2019-01-13 08:05:32
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/uscalfe/uscalfe.h>
#include "mesh_function_utils.h"

namespace lf::uscalfe::test {

struct X {
  X() = default;
  X(int x) : x_(x) {}

  X operator-() const { return X(-x_); }
  bool operator==(const X& rhs) const { return x_ == rhs.x_; }

 private:
  int x_;
};

auto mfScalar = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) { return x[0] * x[0] + x[1]; });

auto mfVector = MeshFunctionGlobal([](const Eigen::Vector2d& x) { return x; });

auto mfVector_dynamic = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) { return Eigen::VectorXd(x); });

auto mfRowVector =
    MeshFunctionGlobal([](auto x) { return Eigen::RowVector2d(x[1], -x[0]); });

auto mfRowVector_dynamic = MeshFunctionGlobal(
    [](auto x) { return Eigen::MatrixXd(Eigen::RowVector2d(x[1], -x[0])); });

auto mfMatrix =
    MeshFunctionGlobal([](auto x) { return (x * x.transpose()).eval(); });

auto mfX = MeshFunctionConstant(X(2));

TEST(meshFunctionUnary, minus) {
  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();

  auto minusScalar = -mfScalar;
  checkMeshFunctionEqual(*mesh, minusScalar, MeshFunctionGlobal([](auto x) {
    return -x[0] * x[0] - x[1];
  }));

  auto minusVector = -mfVector;
  checkMeshFunctionEqual(*mesh, minusVector, MeshFunctionGlobal([](auto x) {
    return (-x).eval();
  }));

  auto minusMatrix = -mfMatrix;
  checkMeshFunctionEqual(*mesh, minusMatrix, MeshFunctionGlobal([](auto x) {
    return (-x * x.transpose()).eval();
  }));

  auto minusX = -mfX;
  checkMeshFunctionEqual(*mesh, minusX, MeshFunctionConstant(X(-2)));
}

TEST(meshFunctionUnary, squaredNorm) {
  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();

  auto snScalar = squaredNorm(mfScalar);
  checkMeshFunctionEqual(*mesh, snScalar, MeshFunctionGlobal([](auto x) {
    return std::pow(x[0] * x[0] + x[1], 2);
  }));

  auto snMfVector = MeshFunctionGlobal([](auto x) { return x.squaredNorm(); });
  checkMeshFunctionEqual(*mesh, squaredNorm(mfVector), snMfVector);
  checkMeshFunctionEqual(*mesh, squaredNorm(mfVector_dynamic), snMfVector);

  auto snMfRowVector =
      MeshFunctionGlobal([](auto x) { return x.squaredNorm(); });
  checkMeshFunctionEqual(*mesh, squaredNorm(mfRowVector), snMfRowVector);
  checkMeshFunctionEqual(*mesh, squaredNorm(mfRowVector_dynamic),
                         snMfRowVector);

  auto snMatrix = squaredNorm(mfMatrix);
  checkMeshFunctionEqual(*mesh, snMatrix, MeshFunctionGlobal([](auto x) {
    return (x * x.transpose()).squaredNorm();
  }));
}

TEST(meshFunctionUnary, transpose) {
  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();

  auto mfVectorT = MeshFunctionGlobal(
      [](auto x) -> Eigen::Matrix<double, 1, 2> { return x.transpose(); });
  checkMeshFunctionEqual(*mesh, transpose(mfVector), mfVectorT);
  checkMeshFunctionEqual(*mesh, transpose(mfVector_dynamic), mfVectorT);

  auto mfRowVectorT = MeshFunctionGlobal(
      [](auto x) -> Eigen::Vector2d { return Eigen::Vector2d(x[1], -x[0]); });
  checkMeshFunctionEqual(*mesh, transpose(mfRowVector), mfRowVectorT);
  checkMeshFunctionEqual(*mesh, transpose(mfRowVector_dynamic), mfRowVectorT);

  checkMeshFunctionEqual(*mesh, transpose(mfMatrix), mfMatrix);
}

}  // namespace lf::uscalfe::test

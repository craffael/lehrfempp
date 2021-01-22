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
#include "mesh_function_utils.h"

namespace lf::mesh::utils::test {

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

auto mfArray = MeshFunctionGlobal(
    [](const Eigen::Array2d& x) -> Eigen::Array2d { return x.array(); });
auto mfArray_dynamic = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) -> Eigen::ArrayXd { return x.array(); });

auto mfRowVector =
    MeshFunctionGlobal([](auto x) { return Eigen::RowVector2d(x[1], -x[0]); });
auto mfRowVector_dynamic = MeshFunctionGlobal(
    [](auto x) { return Eigen::MatrixXd(Eigen::RowVector2d(x[1], -x[0])); });
auto mfRowArray = MeshFunctionGlobal(
    [](auto x) { return Eigen::Array2d(x[1], -x[0]).transpose().eval(); });
auto mfRowArray_dynamic =
    MeshFunctionGlobal([](auto x) -> Eigen::Array<double, 1, Eigen::Dynamic> {
      return Eigen::Array2d(x[1], -x[0]).transpose().eval();
    });

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

  auto minusVRef = MeshFunctionGlobal([](auto x) { return (-x).eval(); });
  checkMeshFunctionEqual(*mesh, -mfVector, minusVRef);
  checkMeshFunctionEqual(*mesh, -mfVector_dynamic, minusVRef);
  checkMeshFunctionEqual(*mesh, -mfArray, minusVRef);
  checkMeshFunctionEqual(*mesh, -mfArray_dynamic, minusVRef);

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
  checkMeshFunctionEqual(*mesh, transpose(mfArray), mfVectorT);
  checkMeshFunctionEqual(*mesh, transpose(mfArray_dynamic), mfVectorT);

  auto mfRowVectorT = MeshFunctionGlobal(
      [](auto x) -> Eigen::Vector2d { return Eigen::Vector2d(x[1], -x[0]); });
  checkMeshFunctionEqual(*mesh, transpose(mfRowVector), mfRowVectorT);
  checkMeshFunctionEqual(*mesh, transpose(mfRowVector_dynamic), mfRowVectorT);
  checkMeshFunctionEqual(*mesh, transpose(mfRowArray), mfRowVectorT);
  checkMeshFunctionEqual(*mesh, transpose(mfRowArray_dynamic), mfRowVectorT);

  checkMeshFunctionEqual(*mesh, transpose(mfMatrix), mfMatrix);
}

TEST(meshFunctionUnary, adjoint) {
  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();
  using c_t = std::complex<double>;
  MeshFunctionConstant a(c_t{1, 1});
  MeshFunctionConstant b(c_t{1, -1});
  auto mfVectorTrans = MeshFunctionGlobal(
      [](auto x) -> Eigen::Matrix<double, 1, 2> { return x.transpose(); });
  checkMeshFunctionEqual(*mesh, adjoint(a * mfVector), b * mfVectorTrans);
  checkMeshFunctionEqual(*mesh, adjoint(a * mfVector_dynamic),
                         b * mfVectorTrans);

  auto mfRowVectorT = MeshFunctionGlobal(
      [](auto x) -> Eigen::Vector2d { return Eigen::Vector2d(x[1], -x[0]); });
  checkMeshFunctionEqual(*mesh, adjoint(a * mfRowVector), b * mfRowVectorT);
  checkMeshFunctionEqual(*mesh, adjoint(a * mfRowVector_dynamic),
                         b * mfRowVectorT);

  checkMeshFunctionEqual(*mesh, adjoint(a * mfMatrix), b * mfMatrix);
}

TEST(meshFunctionUnary, conjugate) {
  auto reader = io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto mesh = reader.mesh();
  using c_t = std::complex<double>;
  MeshFunctionConstant a(c_t{1, 1});
  MeshFunctionConstant b(c_t{1, -1});

  checkMeshFunctionEqual(*mesh, conjugate(mfScalar), mfScalar);
  checkMeshFunctionEqual(*mesh, conjugate(a * mfScalar), b * mfScalar);

  checkMeshFunctionEqual(*mesh, conjugate(a * mfVector), b * mfVector);
  checkMeshFunctionEqual(*mesh, conjugate(a * mfVector_dynamic), b * mfVector);
  checkMeshFunctionEqual(*mesh, conjugate(a * mfArray), b * mfVector);
  checkMeshFunctionEqual(*mesh, conjugate(a * mfArray_dynamic), b * mfVector);

  checkMeshFunctionEqual(*mesh, conjugate(a * mfRowVector), b * mfRowVector);
  checkMeshFunctionEqual(*mesh, conjugate(a * mfRowVector_dynamic),
                         b * mfRowVector);
  checkMeshFunctionEqual(*mesh, conjugate(a * mfRowArray), b * mfRowVector);
  checkMeshFunctionEqual(*mesh, conjugate(a * mfRowArray_dynamic),
                         b * mfRowVector);

  checkMeshFunctionEqual(*mesh, conjugate(a * mfMatrix), b * mfMatrix);
}

}  // namespace lf::mesh::utils::test

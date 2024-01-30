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

#include "lf/mesh/test_utils/test_meshes.h"
#include "mesh_function_utils.h"

namespace lf::mesh::utils {

// Make sure the functionality in MeshFunctionBinary which requires concepts
// (such as MeshFunction), rely indeed only on the functionality in these
// concepts. We do this by instantiating the corresponding functions/classes
// using Archetypes.
template auto operator+(const MeshFunctionAT<double>&,
                        const MeshFunctionAT<double>&)
    -> MeshFunctionBinary<internal::OperatorAddition, MeshFunctionAT<double>,
                          MeshFunctionAT<double>>;

template auto operator-(const MeshFunctionAT<double>&,
                        const MeshFunctionAT<double>&)
    -> MeshFunctionBinary<internal::OperatorSubtraction, MeshFunctionAT<double>,
                          MeshFunctionAT<double>>;

template auto operator*(const MeshFunctionAT<double>&,
                        const MeshFunctionAT<double>&)
    -> MeshFunctionBinary<internal::OperatorMultiplication,
                          MeshFunctionAT<double>, MeshFunctionAT<double>>;

}  // namespace lf::mesh::utils

namespace lf::mesh::utils::test {

class X {
 public:
  X() {}
  X(int x) : x_(x) {}

  X operator+(const X& rhs) const { return X(x_ + rhs.x_); }
  X operator-(const X& rhs) const { return X(x_ - rhs.x_); }
  X operator*(const X& rhs) const { return X(x_ * rhs.x_); };
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
auto mfVectorA_dynamic = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) { return Eigen::VectorXd(x); });
auto mfVectorB_dynamic = MeshFunctionGlobal([](const Eigen::Vector2d& x) {
  return Eigen::VectorXd(Eigen::Vector2d(x[0], 2 * x[0]));
});

auto mfArrayA = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) -> Eigen::Array2d { return x.array(); });
auto mfArrayB =
    MeshFunctionGlobal([](const Eigen::Vector2d& x) -> Eigen::Array2d {
      return Eigen::Array2d(x[0], 2 * x[0]);
    });
auto mfArrayA_dynamic = MeshFunctionGlobal(
    [](const Eigen::Vector2d& x) -> Eigen::ArrayXd { return x.array(); });
auto mfArrayB_dynamic =
    MeshFunctionGlobal([](const Eigen::Vector2d& x) -> Eigen::ArrayXd {
      return Eigen::Array2d(x[0], 2 * x[0]);
    });

auto mfMatrixA =
    MeshFunctionGlobal([](auto x) { return (x * x.transpose()).eval(); });
auto mfMatrixB = MeshFunctionGlobal([](auto x) -> Eigen::Matrix2d {
  return (Eigen::Matrix2d() << 0, 1, x[0], x[1]).finished();
});
auto mfArrayA22 = MeshFunctionGlobal(
    [](auto x) { return (x * x.transpose()).array().eval(); });
auto mfArrayB22 = MeshFunctionGlobal(
    [](auto x) { return (Eigen::Array22d() << 0, 1, x[0], x[1]).finished(); });

auto mfXA = MeshFunctionConstant(X(1));
auto mfXB = MeshFunctionConstant(X(2));

TEST(meshFunctionBinary, Addition) {
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  auto sum = mfA + mfB;
  checkMeshFunctionEqual(*mesh, sum, MeshFunctionGlobal([](auto x) {
    return x[0] * x[0] + 2 * x[1] + x[0];
  }));

  auto mfVecAdd = MeshFunctionGlobal(
      [](auto x) { return Eigen::Vector2d(2 * x[0], x[1] + 2 * x[0]); });
  checkMeshFunctionEqual(*mesh, mfVectorA + mfVectorB, mfVecAdd);
  checkMeshFunctionEqual(*mesh, mfVectorA_dynamic + mfVectorB, mfVecAdd);
  checkMeshFunctionEqual(*mesh, mfVectorA + mfVectorB_dynamic, mfVecAdd);
  checkMeshFunctionEqual(*mesh, mfVectorA_dynamic + mfVectorB_dynamic,
                         mfVecAdd);
  checkMeshFunctionEqual(*mesh, mfArrayA + mfArrayB, mfVecAdd);
  checkMeshFunctionEqual(*mesh, mfArrayA_dynamic + mfArrayB, mfVecAdd);
  checkMeshFunctionEqual(*mesh, mfArrayA + mfArrayB_dynamic, mfVecAdd);
  checkMeshFunctionEqual(*mesh, mfArrayA_dynamic + mfArrayB_dynamic, mfVecAdd);

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
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  auto sub = mfA - mfB;
  checkMeshFunctionEqual(*mesh, sub, MeshFunctionGlobal([](auto x) {
    return x[0] * x[0] - x[0];
  }));

  auto mfVecDif = MeshFunctionGlobal(
      [](auto x) { return Eigen::Vector2d(0., x[1] - 2 * x[0]); });
  checkMeshFunctionEqual(*mesh, mfVectorA - mfVectorB, mfVecDif);
  checkMeshFunctionEqual(*mesh, mfVectorA_dynamic - mfVectorB, mfVecDif);
  checkMeshFunctionEqual(*mesh, mfVectorA - mfVectorB_dynamic, mfVecDif);
  checkMeshFunctionEqual(*mesh, mfVectorA_dynamic - mfVectorB_dynamic,
                         mfVecDif);
  checkMeshFunctionEqual(*mesh, mfArrayA - mfArrayB, mfVecDif);
  checkMeshFunctionEqual(*mesh, mfArrayA_dynamic - mfArrayB, mfVecDif);
  checkMeshFunctionEqual(*mesh, mfArrayA - mfArrayB_dynamic, mfVecDif);
  checkMeshFunctionEqual(*mesh, mfArrayA_dynamic - mfArrayB_dynamic, mfVecDif);

  auto subMatrix = mfMatrixA - mfMatrixB;
  checkMeshFunctionEqual(*mesh, subMatrix, MeshFunctionGlobal([](auto x) {
    return (Eigen::Matrix2d() << x[0] * x[0], x[0] * x[1] - 1,
            x[0] * x[1] - x[0], x[1] * x[1] - x[1])
        .finished();
  }));

  auto xSub = mfXA - mfXB;
  checkMeshFunctionEqual(*mesh, xSub, MeshFunctionConstant(X(-1)));
}

TEST(meshFunctionBinary, Multiplication) {
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  auto mult = mfA * mfB;
  checkMeshFunctionEqual(*mesh, mult, MeshFunctionGlobal([](auto x) {
    return (x[0] * x[0] + x[1]) * (x[0] + x[1]);
  }));

  // mfA * mfVectorA
  auto mfVecMult = MeshFunctionGlobal([](auto x) {
    return Eigen::Vector2d(x[0] * (x[0] * x[0] + x[1]),
                           x[1] * (x[0] * x[0] + x[1]));
  });
  checkMeshFunctionEqual(*mesh, mfA * mfVectorA, mfVecMult);
  checkMeshFunctionEqual(*mesh, mfVectorA * mfA, mfVecMult);
  checkMeshFunctionEqual(*mesh, mfA * mfVectorA_dynamic, mfVecMult);
  checkMeshFunctionEqual(*mesh, mfVectorA_dynamic * mfA, mfVecMult);

  // mfA * mfArrayA
  checkMeshFunctionEqual(*mesh, mfA * mfArrayA, mfVecMult);
  checkMeshFunctionEqual(*mesh, mfArrayA * mfA, mfVecMult);
  checkMeshFunctionEqual(*mesh, mfA * mfArrayA_dynamic, mfVecMult);
  checkMeshFunctionEqual(*mesh, mfArrayA_dynamic * mfA, mfVecMult);

  // transpose(mfVectorA)*mfVectorB
  auto mfVecMult2 = MeshFunctionGlobal([](auto x) {
    return (Eigen::Matrix<double, 1, 1>() << x[0] * x[0] + 2 * x[1] * x[0])
        .finished();
  });
  checkMeshFunctionEqual(*mesh, transpose(mfVectorA) * mfVectorB, mfVecMult2);
  checkMeshFunctionEqual(*mesh, transpose(mfVectorA_dynamic) * mfVectorB,
                         mfVecMult2);
  checkMeshFunctionEqual(*mesh, transpose(mfVectorA) * mfVectorB_dynamic,
                         mfVecMult2);
  checkMeshFunctionEqual(
      *mesh, transpose(mfVectorA_dynamic) * mfVectorB_dynamic, mfVecMult2);
  checkMeshFunctionEqual(*mesh, transpose(mfVectorB) * mfVectorA, mfVecMult2);
  checkMeshFunctionEqual(*mesh, transpose(mfVectorB_dynamic) * mfVectorA,
                         mfVecMult2);
  checkMeshFunctionEqual(*mesh, transpose(mfVectorB) * mfVectorA_dynamic,
                         mfVecMult2);
  checkMeshFunctionEqual(
      *mesh, transpose(mfVectorB_dynamic) * mfVectorA_dynamic, mfVecMult2);

  // mfArrayA*mfArrayB
  auto mfArrayMult = MeshFunctionGlobal(
      [](auto x) { return Eigen::Array2d(x[0] * x[0], 2 * x[0] * x[1]); });
  checkMeshFunctionEqual(*mesh, mfArrayA * mfArrayB, mfArrayMult);
  checkMeshFunctionEqual(*mesh, mfArrayA_dynamic * mfArrayB, mfArrayMult);
  checkMeshFunctionEqual(*mesh, mfArrayA * mfArrayB_dynamic, mfArrayMult);
  checkMeshFunctionEqual(*mesh, mfArrayA_dynamic * mfArrayB_dynamic,
                         mfArrayMult);
  checkMeshFunctionEqual(*mesh, mfArrayB * mfArrayA, mfArrayMult);
  checkMeshFunctionEqual(*mesh, mfArrayB_dynamic * mfArrayA, mfArrayMult);
  checkMeshFunctionEqual(*mesh, mfArrayB * mfArrayA_dynamic, mfArrayMult);
  checkMeshFunctionEqual(*mesh, mfArrayB_dynamic * mfArrayA_dynamic,
                         mfArrayMult);

  // mfMatrixA * mfVectorA
  auto mfMatrixA_mfVectorA = MeshFunctionGlobal([](auto x) {
    return Eigen::Vector2d(x[0] * x[0] * x[0] + x[0] * x[1] * x[1],
                           x[0] * x[0] * x[1] + x[1] * x[1] * x[1]);
  });
  checkMeshFunctionEqual(*mesh, mfMatrixA * mfVectorA, mfMatrixA_mfVectorA);
  checkMeshFunctionEqual(*mesh, mfMatrixA * mfVectorA_dynamic,
                         mfMatrixA_mfVectorA);

  // mfArrayA22 * mfArrayB22
  auto mfArrayA22_mfArrayA = MeshFunctionGlobal([](auto x) -> Eigen::Array22d {
    return (Eigen::Array22d() << 0, x[0] * x[1], x[0] * x[0] * x[1],
            x[1] * x[1] * x[1])
        .finished();
  });
  checkMeshFunctionEqual(*mesh, mfArrayA22 * mfArrayB22, mfArrayA22_mfArrayA);

  // multiplication of two arbitrary types supporting operator*
  auto xSub = mfXA * mfXB;
  checkMeshFunctionEqual(*mesh, mfXA * mfXB, MeshFunctionConstant(X(2)));
}

}  // namespace lf::mesh::utils::test

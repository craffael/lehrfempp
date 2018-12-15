/**
 * @file
 * @brief Test if the lf::fe::isMeshFunction works as expected.
 * @author Raffael Casagrande
 * @date   2018-12-15 04:12:13
 * @copyright MIT License
 */

#include <lf/fe/is_mesh_function.h>

#include <gtest/gtest.h>

namespace lf::fe::test {

struct ProperMeshFunction {
  auto operator()(const mesh::Entity& e, const Eigen::MatrixXd& local) const {
    return std::vector<double>{1.0};
  }
};

struct NotCopyableMeshFunction {
  NotCopyableMeshFunction(const NotCopyableMeshFunction&) = delete;

  auto operator()(const mesh::Entity& e, const Eigen::MatrixXd& local) const {
    return std::vector<double>{1.0};
  }
};

struct NotMoveableMeshFunction {
  NotMoveableMeshFunction(const NotMoveableMeshFunction&) = delete;

  auto operator()(const mesh::Entity& e, const Eigen::MatrixXd& local) const {
    return std::vector<double>{1.0};
  }
};

struct NotCallableMeshFunction {};

TEST(isMeshFunctionTestSuite, basicTests) {
  EXPECT_TRUE(isMeshFunction<ProperMeshFunction>);
  EXPECT_TRUE((isMeshFunction<ProperMeshFunction, double>));
  EXPECT_FALSE((isMeshFunction<ProperMeshFunction, float>));
  EXPECT_TRUE(isMeshFunction<const ProperMeshFunction>);
  EXPECT_TRUE((isMeshFunction<const ProperMeshFunction, double>));
  EXPECT_FALSE((isMeshFunction<const ProperMeshFunction, std::vector<double>>));

  // references are not MeshFunctions
  EXPECT_FALSE(isMeshFunction<ProperMeshFunction&>);
  EXPECT_FALSE((isMeshFunction<ProperMeshFunction&, double>));
  EXPECT_FALSE(isMeshFunction<const ProperMeshFunction&>);
  EXPECT_FALSE((isMeshFunction<const ProperMeshFunction&, double>));
  EXPECT_FALSE(isMeshFunction<ProperMeshFunction&&>);
  EXPECT_FALSE((isMeshFunction<ProperMeshFunction&&, double>));

  // also pointers don't work
  EXPECT_FALSE(isMeshFunction<ProperMeshFunction*>);

  // check with lambda
  auto lambda = [](const mesh::Entity& e, const Eigen::MatrixXd& local) {
    return std::vector<double>{1.0};
  };
  EXPECT_TRUE(isMeshFunction<decltype(lambda)>);
  EXPECT_TRUE((isMeshFunction<decltype(lambda), double>));

  // test that any of the other expressions don't work
  EXPECT_FALSE(isMeshFunction<NotCopyableMeshFunction>);
  EXPECT_FALSE((isMeshFunction<NotCopyableMeshFunction, double>));
  EXPECT_FALSE(isMeshFunction<NotMoveableMeshFunction>);
  EXPECT_FALSE((isMeshFunction<NotMoveableMeshFunction, double>));
  EXPECT_FALSE(isMeshFunction<NotCallableMeshFunction>);
  EXPECT_FALSE((isMeshFunction<NotCallableMeshFunction, double>));

  // if operator() expects wrong number of arguments
  auto l2 = [](const mesh::Entity& e) { return std::vector{1.0}; };
  EXPECT_FALSE(isMeshFunction<decltype(l2)>);
  EXPECT_FALSE((isMeshFunction<decltype(l2), double>));

  // if operator() has a template argument
  auto l3 = [](auto& e, auto local) { return std::vector{1.0}; };
  EXPECT_TRUE(isMeshFunction<decltype(l3)>);
  EXPECT_TRUE((isMeshFunction<decltype(l3), double>));
  EXPECT_FALSE((isMeshFunction<decltype(l3), float>));

  // if operator() expects std::vector instead of Eigen::MatrixXd
  auto l4 = [](const mesh::Entity&, std::vector<double> local) {
    return std::vector{1.0};
  };
  EXPECT_FALSE(isMeshFunction<decltype(l4)>);
  EXPECT_FALSE((isMeshFunction<decltype(l4), double>));

  // if operator() doesn't return std::vector
  auto l5 = [](const mesh::Entity&, Eigen::MatrixXd) { return 1; };
  EXPECT_FALSE(isMeshFunction<decltype(l5)>);
  EXPECT_FALSE((isMeshFunction<decltype(l5), int>));

  // if operator() returns a vector of 3x1 eigen vectors
  auto l6 = [](const mesh::Entity&, const Eigen::MatrixXd&) {
    return std::vector{Eigen::Vector3d(1.0, 1.0, 1.0)};
  };
  EXPECT_TRUE(isMeshFunction<decltype(l6)>);
  EXPECT_TRUE((isMeshFunction<decltype(l6), Eigen::Vector3d>));

  // if operator() accepts Eigen matrix by value
  auto l7 = [](const mesh::Entity&, Eigen::MatrixXd) {
    return std::vector{1.0};
  };
  EXPECT_TRUE(isMeshFunction<decltype(l7)>);
  EXPECT_TRUE((isMeshFunction<decltype(l7), double>));

  // if operator() accepts Eigen matrix by mutable reference
  auto l8 = [](const mesh::Entity&, Eigen::MatrixXd&) {
    return std::vector{1.0};
  };
  EXPECT_FALSE(isMeshFunction<decltype(l8)>);
  EXPECT_FALSE((isMeshFunction<decltype(l8), double>));
}

}  // namespace lf::fe::test

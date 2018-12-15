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
  Eigen::MatrixXd operator()(const mesh::Entity& e,
                             const Eigen::MatrixXd& local) const {
    return Eigen::MatrixXd(3, 3);
  }
};

struct NotCopyableMeshFunction {
  NotCopyableMeshFunction(const NotCopyableMeshFunction&) = delete;

  Eigen::MatrixXd operator()(const mesh::Entity& e,
                             const Eigen::MatrixXd& local) const {
    return Eigen::MatrixXd(3, 3);
  }
};

struct NotMoveableMeshFunction {
  NotMoveableMeshFunction(const NotMoveableMeshFunction&) = delete;

  Eigen::MatrixXd operator()(const mesh::Entity& e,
                             const Eigen::MatrixXd& local) const {
    return Eigen::MatrixXd(3, 3);
  }
};

struct NotCallableMeshFunction {};

TEST(isMeshFunctionTestSuite, basicTests) {
  EXPECT_EQ(lf::fe::isMeshFunction<ProperMeshFunction>, true);
  EXPECT_EQ(lf::fe::isMeshFunction<const ProperMeshFunction>, true);

  // references are not MeshFunctions
  EXPECT_EQ(lf::fe::isMeshFunction<ProperMeshFunction&>, false);
  EXPECT_EQ(lf::fe::isMeshFunction<const ProperMeshFunction&>, false);
  EXPECT_EQ(lf::fe::isMeshFunction<ProperMeshFunction&&>, false);

  // also pointers don't work
  EXPECT_EQ(lf::fe::isMeshFunction<ProperMeshFunction*>, false);

  // check with lambda
  auto lambda = [](const mesh::Entity& e, const Eigen::MatrixXd& local) {
    return Eigen::MatrixXcd(3, 3);
  };
  EXPECT_EQ(lf::fe::isMeshFunction<decltype(lambda)>, true);

  // test if any of the other expressions don't work
  EXPECT_EQ(lf::fe::isMeshFunction<NotCopyableMeshFunction>, false);
  EXPECT_EQ(lf::fe::isMeshFunction<NotMoveableMeshFunction>, false);
  EXPECT_EQ(lf::fe::isMeshFunction<NotCallableMeshFunction>, false);

  // if operator() expects wrong number of arguments
  auto l2 = [](const mesh::Entity& e) { return Eigen::MatrixXd(3, 3); };
  EXPECT_EQ(lf::fe::isMeshFunction<decltype(l2)>, false);

  // if operator() has a template argument
  auto l3 = [](auto& e, auto local) { return Eigen::MatrixXd(3, 3); };
  EXPECT_EQ(lf::fe::isMeshFunction<decltype(l3)>, true);

  // if operator() expect std::vector instead of Eigen::MatrixXd
  auto l4 = [](const mesh::Entity&, std::vector<double> local) {
    return Eigen::MatrixXd(3, 3);
  };
  EXPECT_EQ(lf::fe::isMeshFunction<decltype(l4)>, false);

  // if operator() doesn't return eigen matrix
  auto l5 = [](const mesh::Entity&, Eigen::MatrixXd) {
    return std::vector<double>{1.0};
  };
  EXPECT_FALSE(lf::fe::isMeshFunction<decltype(l5)>);

  // if operator() returns an eigen vector
  auto l6 = [](const mesh::Entity&, const Eigen::MatrixXd&) {
    return Eigen::VectorXd(5);
  };
  EXPECT_TRUE(lf::fe::isMeshFunction<decltype(l6)>);

  // if operator() accepts Eigen matrix by value
  auto l7 = [](const mesh::Entity&, Eigen::MatrixXd) {
    return Eigen::MatrixXd(3, 3);
  };
  EXPECT_TRUE(lf::fe::isMeshFunction<decltype(l7)>);

  // if operator() accepts Eigen matrix by mutable reference
  auto l8 = [](const mesh::Entity&, Eigen::MatrixXd&) {
    return Eigen::MatrixXd(3, 3);
  };
  EXPECT_FALSE(lf::fe::isMeshFunction<decltype(l8)>);
}

}  // namespace lf::fe::test

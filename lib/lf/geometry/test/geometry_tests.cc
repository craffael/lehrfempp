/**
 * @file
 * @brief tests for geometry objects
 * @author Anian Ruoss
 * @date   2018-10-27 15:57:17
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>

using GeometryObjects =
    testing::Types<lf::geometry::TriaO1, lf::geometry::QuadO1>;

static const Eigen::Matrix2d map = (Eigen::Matrix2d() << 1, 2, 3, 4).finished();
static const Eigen::Vector2d offset = (Eigen::Vector2d() << 1, 2).finished();

static Eigen::Vector2d affineMap(const Eigen::Vector2d &x) {
  return map * x + offset;
}

template <class T>
class GeometryTest : public testing::Test {
 protected:
  std::shared_ptr<T> object_;
  std::shared_ptr<Eigen::MatrixXd> refElCoords_;

  GeometryTest() {
    std::shared_ptr<lf::base::RefEl> refEl;

    if (std::is_same<lf::geometry::Point, T>::value) {
      refEl = std::make_shared<lf::base::RefEl>(lf::base::RefEl::kPoint());
    } else if (std::is_same<lf::geometry::SegmentO1, T>::value) {
      refEl = std::make_shared<lf::base::RefEl>(lf::base::RefEl::kSegment());
    } else if (std::is_same<lf::geometry::TriaO1, T>::value) {
      refEl = std::make_shared<lf::base::RefEl>(lf::base::RefEl::kTria());
    } else if (std::is_same<lf::geometry::QuadO1, T>::value) {
      refEl = std::make_shared<lf::base::RefEl>(lf::base::RefEl::kQuad());
    }

    refElCoords_ = std::make_shared<Eigen::MatrixXd>(refEl->NodeCoords());
    Eigen::MatrixXd coords(2, refElCoords_->cols());

    for (int j = 0; j < refElCoords_->cols(); ++j) {
      coords.col(j) = affineMap(refElCoords_->col(j));
    }

    object_ = std::make_shared<T>(coords);
  }
};

TYPED_TEST_CASE(GeometryTest, GeometryObjects);

TYPED_TEST(GeometryTest, checkJacobian) {
  Eigen::MatrixXd jacobian = this->object_->Jacobian(*this->refElCoords_);

  for (int i = 0; i < this->refElCoords_->cols(); ++i) {
    EXPECT_EQ(jacobian.block(0, i * 2, 2, 2), map)
        << "Jacobian incorrect at vertex " << i;
  }
}

TYPED_TEST(GeometryTest, checkJacobianInverseGramian) {
  Eigen::MatrixXd jacInvGram =
      this->object_->JacobianInverseGramian(*this->refElCoords_);

  for (int i = 0; i < this->refElCoords_->cols(); ++i) {
    EXPECT_EQ(jacInvGram.block(0, i * 2, 2, 2),
              map * (map.transpose() * map).inverse())
        << "JacobianInverseGramian incorrect at vertex " << i;
  }
}

TYPED_TEST(GeometryTest, checkIntegrationElement) {
  Eigen::VectorXd integrationElement =
      this->object_->IntegrationElement(*this->refElCoords_);

  for (int i = 0; i < this->refElCoords_->cols(); ++i) {
    EXPECT_EQ(integrationElement(i),
              std::sqrt((map.transpose() * map).determinant()))
        << "IntegrationElement incorrect at vertex " << i;
  }
}

/**
 * @file
 * @brief Test BrepMeshFactoryTransfinite
 * @author Raffael Casagrande
 * @date   2021-02-11 02:00:28
 * @copyright MIT License
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <lf/brep/geom/geom.h>
#include <lf/brep/test_utils/curve_circle.h>
#include <lf/brep/test_utils/fake_brep_model.h>
#include <lf/quad/quad.h>

namespace lf::brep::geom::test {

class MockMeshFactory : public mesh::MeshFactory {
  int point_count = 0;

 public:
  MOCK_METHOD(dim_t, DimWorld, (), (const, override));
  MOCK_METHOD(dim_t, DimMesh, (), (const, override));
  MOCK_METHOD(size_type, AddPoint, (Eigen::VectorXd coord), (override));
  MOCK_METHOD(size_type, AddPoint,
              (std::unique_ptr<geometry::Geometry> && geometry), (override));
  MOCK_METHOD(size_type, AddEntity,
              (base::RefEl ref_el, const nonstd::span<const size_type>& nodes,
               std::unique_ptr<geometry::Geometry>&& geometry),
              (override));
  MOCK_METHOD(std::shared_ptr<mesh::Mesh>, Build, (), (override));
};

using ::testing::_;
using ::testing::An;
using ::testing::Return;
using ::testing::SaveArg;

TEST(lf_brep_geom, BrepMeshFactoryTransfiniteSegmentTest) {
  auto brep_model = std::make_shared<test_utils::FakeBrepModel>();
  brep_model->AddCurve(
      std::make_shared<test_utils::CurveCircle>(Eigen::Vector3d(0, 0, 0), 1.0));
  brep_model->AddCurve(std::make_shared<test_utils::CurveCircle>(
      Eigen::Vector3d(-1, -1, 0), 1.));

  auto mmf_ = std::make_unique<MockMeshFactory>();
  auto& mmf = *mmf_.get();

  EXPECT_CALL(mmf, DimWorld).WillRepeatedly(Return(2));
  EXPECT_CALL(mmf, DimMesh).WillRepeatedly(Return(2));

  BrepMeshFactoryTransfinite mf(std::move(mmf_), brep_model);

  // add seven points:
  int point_count = 0;
  EXPECT_CALL(mmf, AddPoint(An<Eigen::VectorXd>()))
      .Times(7)
      .WillRepeatedly([&](auto) { return point_count++; });
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 0)), 0);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(1, 0)), 1);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 1)), 2);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 0.5)), 3);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(-1, -2)), 4);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(-1 / std::sqrt(2), 1 / std::sqrt(2))),
            5);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(-1 / std::sqrt(2), -1 / std::sqrt(2))),
            6);

  // add a segment which doesn't lie on the boundary:
  auto segment = std::make_unique<geometry::SegmentO1>(
      (Eigen::MatrixXd(2, 2) << 0, 0, 0, 0.5).finished());
  auto segment_ptr = segment.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 0);
        EXPECT_EQ(nodes[1], 3);
        EXPECT_EQ(nodes.size(), 2);
        EXPECT_EQ(a.get(), segment_ptr);
        return 0;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kSegment(),
                   std::array<base::size_type, 2>{0, 3}, std::move(segment)),
      0);

  // add a segment which has one point on the boundary:
  segment = std::make_unique<geometry::SegmentO1>(
      (Eigen::MatrixXd(2, 2) << 0, 1, 0, 0).finished());
  segment_ptr = segment.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 0);
        EXPECT_EQ(nodes[1], 1);
        EXPECT_EQ(nodes.size(), 2);
        EXPECT_EQ(a.get(), segment_ptr);
        return 1;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kSegment(),
                   std::array<base::size_type, 2>{0, 1}, std::move(segment)),
      1);

  // add a segment that has both points on the boundary:
  segment = std::make_unique<geometry::SegmentO1>(
      (Eigen::MatrixXd(2, 2) << 1, 0, 0, 1).finished());
  segment_ptr = segment.get();
  auto qr = quad::make_QuadRule(base::RefEl::kSegment(), 20);
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 1);
        EXPECT_EQ(nodes[1], 2);
        EXPECT_EQ(nodes.size(), 2);
        EXPECT_NE(a.get(), segment_ptr);
        EXPECT_EQ(a->DimGlobal(), 2);

        auto global = a->Global(qr.Points());
        for (int i = 0; i < global.cols(); ++i) {
          EXPECT_NEAR(global.col(i).norm(), 1., 1e-6);
        }
        EXPECT_NEAR(a->IntegrationElement(qr.Points()).dot(qr.Weights()),
                    base::kPi / 2, 1e-6);

        return 2;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kSegment(),
                   std::array<base::size_type, 2>{1, 2}, std::move(segment)),
      2);

  // lets  add a segment that has a point on each circle:
  segment = std::make_unique<geometry::SegmentO1>(
      (Eigen::MatrixXd(2, 2) << -1, 1, -2, 0).finished());
  segment_ptr = segment.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 4);
        EXPECT_EQ(nodes[1], 1);
        EXPECT_EQ(nodes.size(), 2);
        EXPECT_EQ(a.get(), segment_ptr);
        return 3;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kSegment(),
                   std::array<base::size_type, 2>{4, 1}, std::move(segment)),
      3);

  // lets add a segment that goes over the periodic boundary:
  double sqrt2 = std::sqrt(2.);
  segment = std::make_unique<geometry::SegmentO1>(
      (Eigen::MatrixXd(2, 2) << -1 / sqrt2, -1 / sqrt2, 1 / sqrt2, -1 / sqrt2)
          .finished());
  segment_ptr = segment.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 5);
        EXPECT_EQ(nodes[1], 6);
        EXPECT_EQ(nodes.size(), 2);
        EXPECT_NE(a.get(), segment_ptr);
        EXPECT_EQ(a->DimGlobal(), 2);

        auto global = a->Global(qr.Points());
        for (int i = 0; i < global.cols(); ++i) {
          EXPECT_NEAR(global.col(i).norm(), 1., 1e-6);
        }
        EXPECT_NEAR(a->IntegrationElement(qr.Points()).dot(qr.Weights()),
                    base::kPi / 2, 1e-6);

        return 4;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kSegment(),
                   std::array<base::size_type, 2>{5, 6}, std::move(segment)),
      4);
}

TEST(lf_brep_geom, BrepMeshFactoryTransfiniteTriaDeathTest) {
  auto brep_model = std::make_shared<test_utils::FakeBrepModel>();
  brep_model->AddCurve(
      std::make_shared<test_utils::CurveCircle>(Eigen::Vector3d(0, 0, 0), 1.0));
  auto origin2 = Eigen::Vector3d(0, -1, 0.);
  brep_model->AddCurve(std::make_shared<test_utils::CurveCircle>(origin2, 2));

  auto mmf_ = std::make_unique<MockMeshFactory>();
  auto& mmf = *mmf_.get();

  EXPECT_CALL(mmf, DimWorld).WillRepeatedly(Return(2));
  EXPECT_CALL(mmf, DimMesh).WillRepeatedly(Return(2));

  BrepMeshFactoryTransfinite mf(std::move(mmf_), brep_model);

  // add six points:
  int point_count = 0;
  double sqrt2 = std::sqrt(2);
  EXPECT_CALL(mmf, AddPoint(An<Eigen::VectorXd>()))
      .Times(9)
      .WillRepeatedly([&](auto) { return point_count++; });
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 0)), 0);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(1, 0)), 1);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 1)), 2);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 0.5)), 3);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0.5, 0)), 4);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, -1)), 5);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(-1 / sqrt2, 1 / sqrt2)), 6);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(-1 / sqrt(2), -1 / sqrt(2))), 7);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(-2, -1)), 8);

  // add a tria which doesn't lie on the boundary:
  auto tria = std::make_unique<geometry::TriaO1>(
      (Eigen::MatrixXd(2, 3) << 0, 0.5, 0, 0, 0, 0.5).finished());
  auto* tria_ptr = tria.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 0);
        EXPECT_EQ(nodes[1], 4);
        EXPECT_EQ(nodes[2], 3);
        EXPECT_EQ(nodes.size(), 3);
        EXPECT_EQ(a.get(), tria_ptr);
        return 0;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kTria(),
                   std::array<base::size_type, 3>{0, 4, 3}, std::move(tria)),
      0);

  // add a tria which has one node on the boundary:
  tria = std::make_unique<geometry::TriaO1>(
      (Eigen::MatrixXd(2, 3) << 0, 1., 0, 0, 0, 0.5).finished());
  tria_ptr = tria.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 0);
        EXPECT_EQ(nodes[1], 1);
        EXPECT_EQ(nodes[2], 3);
        EXPECT_EQ(nodes.size(), 3);
        EXPECT_EQ(a.get(), tria_ptr);
        return 1;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kTria(),
                   std::array<base::size_type, 3>{0, 1, 3}, std::move(tria)),
      1);

  // add a tria which has two nodes on the boundary:
  tria = std::make_unique<geometry::TriaO1>(
      (Eigen::MatrixXd(2, 3) << 0, 1., 0, 0, 0, 1.).finished());
  tria_ptr = tria.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 0);
        EXPECT_EQ(nodes[1], 1);
        EXPECT_EQ(nodes[2], 2);
        EXPECT_EQ(nodes.size(), 3);
        EXPECT_NE(a.get(), tria_ptr);

        auto global = a->SubGeometry(1, 0)->Global(
            Eigen::RowVectorXd::LinSpaced(11, 0, 1));
        EXPECT_TRUE(global.row(1).norm() < 1e-6);

        global = a->SubGeometry(1, 1)->Global(
            Eigen::RowVectorXd::LinSpaced(11, 0, 1));
        EXPECT_TRUE(
            global.colwise().norm().isApprox(Eigen::RowVectorXd::Ones(11)));

        global = a->SubGeometry(1, 2)->Global(
            Eigen::RowVectorXd::LinSpaced(11, 0, 1));
        EXPECT_TRUE(global.row(0).norm() < 1e-6);

        return 2;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kTria(),
                   std::array<base::size_type, 3>{0, 1, 2}, std::move(tria)),
      2);

  // add tria which has three nodes on the same boundary
  // -> should throw an error because then the tria has volume 0
  tria = std::make_unique<geometry::TriaO1>(
      (Eigen::MatrixXd(2, 3) << 0, 1., 0, -1, 0, 1.).finished());
  tria_ptr = tria.get();
  ASSERT_DEATH(
      mf.AddEntity(base::RefEl::kTria(),
                   std::array<base::size_type, 3>{5, 1, 2}, std::move(tria)),
      "volume 0");

  // add a tria where one segment is on circle 1 and one on circle 0
  tria = std::make_unique<geometry::TriaO1>(
      (Eigen::MatrixXd(2, 3) << 1, 0., -2, 0, 1, -1).finished());
  tria_ptr = tria.get();
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes.size(), 3);
        EXPECT_NE(a.get(), tria_ptr);

        auto global = a->SubGeometry(1, 0)->Global(
            Eigen::RowVectorXd::LinSpaced(11, 0, 1));
        EXPECT_TRUE(
            global.colwise().norm().isApprox(Eigen::RowVectorXd::Ones(11)));

        global = a->SubGeometry(1, 1)->Global(
            Eigen::RowVectorXd::LinSpaced(11, 0, 1));
        EXPECT_TRUE((global - origin2.topRows(2).replicate(1, 11))
                        .colwise()
                        .norm()
                        .isApprox(Eigen::RowVectorXd::Constant(11, 2)));

        global = a->SubGeometry(1, 2)->Global(
            Eigen::RowVectorXd::LinSpaced(11, 0, 1));
        EXPECT_TRUE(
            global.row(0).isApprox(Eigen::RowVectorXd::LinSpaced(11, -2, 1)));
        EXPECT_TRUE(
            global.row(1).isApprox(Eigen::RowVectorXd::LinSpaced(11, -1, 0)));

        return 3;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kTria(),
                   std::array<base::size_type, 3>{1, 2, 8}, std::move(tria)),
      3);

  // add a tria where one segment goes over the periodic boundary:
  tria = std::make_unique<geometry::TriaO1>(
      (Eigen::MatrixXd(2, 3) << -1 / sqrt2, 0, -1 / sqrt2, -1 / sqrt2, 0,
       1. / sqrt2)
          .finished());
  tria_ptr = tria.get();
  auto qr = quad::make_QuadRule(base::RefEl::kTria(), 10);
  EXPECT_CALL(mmf, AddEntity)
      .WillOnce([&](auto ref_el, auto& nodes,
                    std::unique_ptr<geometry::Geometry>&& a) -> int {
        EXPECT_EQ(nodes[0], 7);
        EXPECT_EQ(nodes[1], 0);
        EXPECT_EQ(nodes[2], 6);
        EXPECT_EQ(nodes.size(), 3);
        EXPECT_NE(a.get(), tria_ptr);

        EXPECT_NEAR(a->IntegrationElement(qr.Points()).dot(qr.Weights()),
                    base::kPi / 4, 1e-6);

        return 4;
      });
  EXPECT_EQ(
      mf.AddEntity(base::RefEl::kTria(),
                   std::array<base::size_type, 3>{7, 0, 6}, std::move(tria)),
      4);
}
}  // namespace lf::brep::geom::test

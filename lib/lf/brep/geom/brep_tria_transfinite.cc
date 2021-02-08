/**
 * @file
 * @brief Implementation of BrepTriaTransfinite
 * @author Raffael Casagrande
 * @date   2021-02-03 04:06:48
 * @copyright MIT License
 */

#include "brep_tria_transfinite.h"

#include <Eigen/src/QR/ColPivHouseholderQR.h>

#include <Eigen/Dense>

#include "brep_curve.h"

namespace lf::brep::geom {

class BrepTriaTransfiniteCurve : public interface::BrepCurve {
 public:
  BrepTriaTransfiniteCurve(const BrepTriaTransfinite& base_tria,
                           const Eigen::Matrix2d& curve_ends)
      : base_tria_(base_tria), curve_ends_(curve_ends) {}

  [[nodiscard]] base::dim_t DimGlobal() const override {
    return base_tria_.DimGlobal();
  }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd GlobalMulti(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(local.rows() == 1, "Expected exactly 1 row.");
    return base_tria_.Global(curve_ends_.col(0) * (1 - local.array()).matrix() +
                             curve_ends_.col(1) * local);
  }
  [[nodiscard]] Eigen::MatrixXd JacobianMulti(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(local.rows() == 1, "Expected exactly 1 row.");
    auto jac =
        base_tria_.Jacobian(curve_ends_.col(0) * (1 - local.array()).matrix() +
                            curve_ends_.col(1) * local);
    Eigen::MatrixXd result(DimGlobal(), local.cols());
    for (int i = 0; i < local.cols(); ++i) {
      result.col(i) = jac.block(0, 2 * i, DimGlobal(), 2) *
                      (curve_ends_.col(1) - curve_ends_.col(0));
    }
    return result;
  }
  [[nodiscard]] std::vector<bool> IsInBoundingBoxMulti(
      const Eigen::MatrixXd& global) const override {
    LF_VERIFY_MSG(false, "not implemented.");
  }
  [[nodiscard]] Eigen::Vector3d GlobalSingle(double local) const override {
    return GlobalMulti((Eigen::Matrix<double, 1, 1>() << local).finished());
  }
  [[nodiscard]] Eigen::Vector3d JacobianSingle(double local) const override {
    return JacobianMulti((Eigen::Matrix<double, 1, 1>() << local).finished());
  }
  [[nodiscard]] std::pair<double, double> Project(
      const Eigen::Vector3d& global) const override {
    LF_VERIFY_MSG(false, "not implemented.");
  }
  [[nodiscard]] bool IsInBoundingBoxSingle(
      const Eigen::Vector3d& global) const override {
    LF_VERIFY_MSG(false, "Not implemented.");
  }
  [[nodiscard]] bool IsInside(double local) const override {
    LF_VERIFY_MSG(false, "not implemented.");
  }

 private:
  BrepTriaTransfinite base_tria_;
  Eigen::Matrix2d curve_ends_;
};

BrepTriaTransfinite::BrepTriaTransfinite(
    std::array<std::pair<const interface::BrepCurve*, Eigen::RowVector2d>, 3>
        curves,
    std::array<bool, 3> delete_curves)
    : curves_(curves), delete_curves_(delete_curves) {
  LF_ASSERT_MSG(
      curves_[0].first->DimGlobal() == curves_[1].first->DimGlobal() &&
          curves_[0].first->DimGlobal() == curves_[2].first->DimGlobal(),
      "Inconsistent global dimension.");
  // Make sure the curves agree at the endpoints:
  LF_ASSERT_MSG(
      curves_[0]
          .first->GlobalSingle(curves_[0].second(1))
          .isApprox(curves_[1].first->GlobalSingle(curves_[1].second(0))),
      "curves don't match in node 1: \n"
          << curves_[0].first->GlobalSingle(curves_[0].second(1)).transpose()
          << " <-> "
          << curves_[1].first->GlobalSingle(curves_[1].second(0)).transpose());
  LF_ASSERT_MSG(
      curves_[1]
          .first->GlobalSingle(curves_[1].second(1))
          .isApprox(curves_[2].first->GlobalSingle(curves_[2].second(0))),
      "curves don't match in node 2:\n"
          << curves_[1].first->GlobalSingle(curves_[1].second(1)).transpose()
          << " <-> "
          << curves_[2].first->GlobalSingle(curves_[2].second(0)).transpose());
  LF_ASSERT_MSG(
      curves_[2]
          .first->GlobalSingle(curves_[2].second(1))
          .isApprox(curves_[0].first->GlobalSingle(curves_[0].second(0))),
      "curves don't match in node 0:\n"
          << curves_[2].first->GlobalSingle(curves_[2].second(1)).transpose()
          << " <-> "
          << curves_[0].first->GlobalSingle(curves_[0].second(0)).transpose());
  node0_ = curves_[0].first->GlobalSingle(curves_[0].second(0));
}

Eigen::MatrixXd BrepTriaTransfinite::Global(
    const Eigen::MatrixXd& local) const {
  auto g0 = [&](auto& x) {
    return curves_[0]
        .first
        ->GlobalMulti((curves_[0].second(0) +
                       (curves_[0].second(1) - curves_[0].second(0)) * x)
                          .matrix())
        .array()
        .eval();
  };
  auto g1 = [&](auto& x) {
    return curves_[1]
        .first
        ->GlobalMulti((curves_[1].second(0) +
                       (curves_[1].second(1) - curves_[1].second(0)) * x)
                          .matrix())
        .array()
        .eval();
  };
  auto g2 = [&](auto& x) {
    return curves_[2]
        .first
        ->GlobalMulti((curves_[2].second(0) +
                       (curves_[2].second(1) - curves_[2].second(0)) * x)
                          .matrix())
        .array()
        .eval();
  };
  auto x = local.row(0).array().eval();
  auto y = local.row(1).array().eval();
  auto xr = x.replicate(DimGlobal(), 1).eval();
  auto yr = y.replicate(DimGlobal(), 1).eval();

  auto t1 = (g0(x) * (1 - yr) - xr * g0(1 - y)).eval();
  auto t2 =
      ((node0_ * (x + y - 1).matrix()).array() + g2(1 - y) * (1 - xr)).eval();
  auto t21 = ((node0_ * (x + y - 1).matrix()).array()).eval();
  auto t22 = (g2(1 - y) * (1 - xr)).eval();
  auto t3 = (yr * g2(x) + xr * g1(y) + yr * g1(1 - x)).eval();

  return (g0(x) * (1 - yr) - xr * g0(1 - y) +
          (node0_ * (x + y - 1).matrix()).array() + g2(1 - y) * (1 - xr) -
          yr * g2(x) + xr * g1(y) + yr * g1(1 - x))
      .matrix();
}

Eigen::MatrixXd BrepTriaTransfinite::Jacobian(
    const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(DimGlobal(), 2 * local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    auto g0 = [&](double x) {
      return curves_[0].first->GlobalSingle(
          curves_[0].second(0) +
          (curves_[0].second(1) - curves_[0].second(0)) * x);
    };
    auto g1 = [&](double x) {
      return curves_[1].first->GlobalSingle(
          curves_[1].second(0) +
          (curves_[1].second(1) - curves_[1].second(0)) * x);
    };
    auto g2 = [&](double x) {
      return curves_[2].first->GlobalSingle(
          curves_[2].second(0) +
          (curves_[2].second(1) - curves_[2].second(0)) * x);
    };
    auto Dg0 = [&](double x) -> Eigen::Vector3d {
      return curves_[0].first->JacobianSingle(
                 curves_[0].second(0) +
                 (curves_[0].second(1) - curves_[0].second(0)) * x) *
             (curves_[0].second(1) - curves_[0].second(0));
    };
    auto Dg1 = [&](double x) -> Eigen::Vector3d {
      return curves_[1].first->JacobianSingle(
                 curves_[1].second(0) +
                 (curves_[1].second(1) - curves_[1].second(0)) * x) *
             (curves_[1].second(1) - curves_[1].second(0));
    };
    auto Dg2 = [&](double x) -> Eigen::Vector3d {
      return curves_[2].first->JacobianSingle(
                 curves_[2].second(0) +
                 (curves_[2].second(1) - curves_[2].second(0)) * x) *
             (curves_[2].second(1) - curves_[2].second(0));
    };
    auto x = local(0, i);
    auto y = local(1, i);

    result.col(2 * i) = node0_ - g0(1 - y) + g1(y) - g2(1 - y) -
                        (y - 1) * Dg0(x) - y * Dg1(1 - x) - y * Dg2(x);
    result.col(2 * i + 1) = node0_ - g0(x) + g1(1 - x) - g2(x) +
                            x * Dg0(1 - y) + x * Dg1(y) + (x - 1) * Dg2(1 - y);
  }
  return result;
}

Eigen::MatrixXd BrepTriaTransfinite::JacobianInverseGramian(
    const Eigen::MatrixXd& local) const {
  auto jac = Jacobian(local);
  Eigen::MatrixXd result(DimGlobal(), 2 * local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    result.block(0, 2 * i, DimGlobal(), 2) =
        (jac.block(0, 2 * i, DimGlobal(), 2).transpose() *
         jac.block(0, 2 * i, DimGlobal(), 2))
            .colPivHouseholderQr()
            .solve(jac.block(0, 2 * i, DimGlobal(), 2).transpose())
            .transpose();
  }
  return result;
}

Eigen::VectorXd BrepTriaTransfinite::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  auto jac = Jacobian(local);
  Eigen::RowVectorXd result(local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    result(i) =
        std::sqrt(std::abs((jac.block(0, 2 * i, DimGlobal(), 2).transpose() *
                            jac.block(0, 2 * i, DimGlobal(), 2))
                               .determinant()));
  }
  return result;
}

std::unique_ptr<geometry::Geometry> BrepTriaTransfinite::SubGeometry(
    dim_t codim, dim_t i) const {
  LF_ASSERT_MSG(codim >= 0 && codim <= 2, "codim out of range.");
  if (codim == 0) {
    LF_ASSERT_MSG(i == 0, "subidx out of range.");
    return std::make_unique<BrepTriaTransfinite>(*this);
  } else if (codim == 1) {
    LF_ASSERT_MSG(i >= 0 && i < 3, "subidx out of range.");
    return std::make_unique<BrepCurve>(curves_[i].first, curves_[i].second,
                                       false);
  }

  // codim==2
  LF_ASSERT_MSG(i >= 0 && i < 3, "subidx out of range.");
  if (i == 0) {
    return std::make_unique<geometry::Point>(node0_);
  } else {
    return std::make_unique<geometry::Point>(
        curves_[i].first->GlobalSingle(curves_[i].second(0)));
  }
}

std::vector<std::unique_ptr<geometry::Geometry>>
BrepTriaTransfinite::ChildGeometry(const geometry::RefinementPattern& ref_pat,
                                   lf::base::dim_t codim) const {
  LF_ASSERT_MSG(codim >= 0 && codim < 3, "codim out of bounds.");
  std::vector<std::unique_ptr<geometry::Geometry>> result;
  result.reserve(ref_pat.NumChildren(codim));
  double lattice_const = static_cast<double>(ref_pat.LatticeConst());
  if (codim == 0) {
    for (auto& c : ref_pat.ChildPolygons(0)) {
      LF_ASSERT_MSG(
          c.cols() == 3,
          "So far a triangle can only be split into other triangles...");
      // construct curves of the triangle:
      Eigen::Matrix2d curve_ends0;
      curve_ends0.col(0) = c.col(0).cast<double>() / lattice_const;
      curve_ends0.col(1) = c.col(1).cast<double>() / lattice_const;
      auto* curve0 = new BrepTriaTransfiniteCurve(*this, curve_ends0);

      Eigen::Matrix2d curve_ends1;
      curve_ends1.col(0) = c.col(1).cast<double>() / lattice_const;
      curve_ends1.col(1) = c.col(2).cast<double>() / lattice_const;
      auto* curve1 = new BrepTriaTransfiniteCurve(*this, curve_ends1);

      Eigen::Matrix2d curve_ends2;
      curve_ends2.col(0) = c.col(2).cast<double>() / lattice_const;
      curve_ends2.col(1) = c.col(0).cast<double>() / lattice_const;
      auto* curve2 = new BrepTriaTransfiniteCurve(*this, curve_ends2);

      Eigen::RowVector2d curve_coords(0., 1.);

      using p_t = std::pair<const interface::BrepCurve*, Eigen::RowVector2d>;
      result.push_back(std::make_unique<BrepTriaTransfinite>(
          std::array<p_t, 3>{p_t{curve0, curve_coords},
                             p_t{curve1, curve_coords},
                             p_t{curve2, curve_coords}},
          std::array<bool, 3>{true, true, true}));
    }
  } else if (codim == 1) {
    for (auto& c : ref_pat.ChildPolygons(1)) {
      LF_ASSERT_MSG(c.cols() == 2, "Expected 2 columns...");
      result.push_back(std::make_unique<BrepCurve>(
          new BrepTriaTransfiniteCurve(*this, c.cast<double>() / lattice_const),
          Eigen::RowVector2d(0., 1.), true));
    }
  } else if (codim == 2) {
    for (auto& c : ref_pat.ChildPolygons(2)) {
      LF_ASSERT_MSG(c.cols() == 2, "Expected 1 column...");
      result.push_back(std::make_unique<geometry::Point>(
          Global(c.cast<double>() / lattice_const)));
    }
  }
  return result;
}

BrepTriaTransfinite::~BrepTriaTransfinite() {
  if (delete_curves_[0]) delete curves_[0].first;
  if (delete_curves_[1]) delete curves_[1].first;
  if (delete_curves_[2]) delete curves_[2].first;
}
}  // namespace lf::brep::geom

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

class BrepTriaTransfiniteCurve : public interface::BrepGeometry {
 public:
  BrepTriaTransfiniteCurve(const BrepTriaTransfinite& base_tria,
                           const Eigen::Matrix2d& curve_ends)
      : base_tria_(base_tria), curve_ends_(curve_ends) {}

  [[nodiscard]] base::dim_t DimGlobal() const override {
    return base_tria_.DimGlobal();
  }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(local.rows() == 1, "Expected exactly 1 row.");
    return base_tria_.Global(curve_ends_.col(0) * (1 - local.array()).matrix() +
                             curve_ends_.col(1) * local);
  }
  [[nodiscard]] Eigen::MatrixXd Jacobian(
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
  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const override {
    LF_VERIFY_MSG(false, "not implemented.");
  }

  [[nodiscard]] std::pair<double, Eigen::VectorXd> Project(
      const Eigen::VectorXd& global) const override {
    LF_VERIFY_MSG(false, "not implemented.");
  }
  [[nodiscard]] bool IsInside(const Eigen::VectorXd& local) const override {
    LF_VERIFY_MSG(false, "not implemented.");
  }

  [[nodiscard]] Eigen::VectorXd Periods() const override {
    LF_VERIFY_MSG(false, "not implemented.");
  }

 private:
  BrepTriaTransfinite base_tria_;
  Eigen::Matrix2d curve_ends_;
};

class SegmentIn2DFather : public geometry::Geometry {
 public:
  SegmentIn2DFather(std::unique_ptr<geometry::Geometry>&& father_geom,
                    const Eigen::Matrix2d& nodes_in_father)
      : father_(std::move(father_geom)) {
    LF_ASSERT_MSG(father_->DimLocal() == 2, "");
    transform_ = nodes_in_father.col(1) - nodes_in_father.col(0);
    offset_ = nodes_in_father.col(0);
  }

  SegmentIn2DFather(const SegmentIn2DFather& other)
      : father_(other.father_->SubGeometry(0, 0)),
        transform_(other.transform_),
        offset_(other.offset_) {}

  [[nodiscard]] dim_t DimLocal() const override { return 1; }
  [[nodiscard]] dim_t DimGlobal() const override {
    return father_->DimGlobal();
  }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kSegment();
  }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override {
    return father_->Global(transform_ * local +
                           offset_.replicate(1, local.cols()));
  }
  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override {
    auto jac = father_->Jacobian(transform_ * local +
                                 offset_.replicate(1, local.cols()));
    Eigen::MatrixXd result(DimGlobal(), local.cols());
    for (int i = 0; i < local.cols(); ++i) {
      result.col(i) = jac.block(0, 2 * i, DimGlobal(), 2) * transform_;
    }
    return result;
  }
  [[nodiscard]] Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override {
    auto jac = Jacobian(local);
    jac.array() /=
        jac.colwise().squaredNorm().array().replicate(DimGlobal(), 1);
    return jac;
  }
  [[nodiscard]] Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override {
    auto jac = Jacobian(local);
    return jac.colwise().norm().transpose();
  }
  [[nodiscard]] std::unique_ptr<Geometry> SubGeometry(dim_t codim,
                                                      dim_t i) const override {
    if (codim == 1) {
      return std::make_unique<geometry::Point>(
          Global(base::RefEl::kSegment().NodeCoords().col(i)));
    }
    LF_ASSERT_MSG(codim == 0, "");
    return std::make_unique<SegmentIn2DFather>(*this);
  }
  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const geometry::RefinementPattern& ref_pat,
      lf::base::dim_t codim) const override {
    std::vector<std::unique_ptr<Geometry>> result;
    result.reserve(ref_pat.NumChildren(codim));
    double lattice_const = ref_pat.LatticeConst();
    if (codim == 1) {
      for (auto& n : ref_pat.ChildPolygons(1)) {
        result.push_back(std::make_unique<geometry::Point>(
            Global(n.cast<double>() / lattice_const)));
      }
    } else {
      for (auto& e : ref_pat.ChildPolygons(0)) {
        LF_ASSERT_MSG(e.cols() == 2, "unexpected number of columns.");
        auto nodes_in_father = (transform_ * e.cast<double>() / lattice_const +
                                offset_.replicate(1, 2))
                                   .eval();
        result.push_back(std::make_unique<SegmentIn2DFather>(
            father_->SubGeometry(0, 0), nodes_in_father));
      }
    }
    return result;
  }

 private:
  std::unique_ptr<geometry::Geometry> father_;
  Eigen::Vector2d transform_;
  Eigen::Vector2d offset_;
};

class TriaInFather : public geometry::Geometry {
 public:
  TriaInFather(std::unique_ptr<Geometry>&& father_geom,
               const Eigen::Matrix<double, 2, 3>& nodes_in_father)
      : father_geom_(std::move(father_geom)) {
    LF_ASSERT_MSG(father_geom_->DimLocal() == 2, "");
    transform_.col(0) = nodes_in_father.col(1) - nodes_in_father.col(0);
    transform_.col(1) = nodes_in_father.col(2) - nodes_in_father.col(0);
    offset_ = nodes_in_father.col(0);
  }

  TriaInFather(const TriaInFather& other)
      : father_geom_(other.father_geom_->SubGeometry(0, 0)),
        transform_(other.transform_),
        offset_(other.offset_) {}

  [[nodiscard]] dim_t DimLocal() const override { return 2; }
  [[nodiscard]] dim_t DimGlobal() const override {
    return father_geom_->DimGlobal();
  }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override {
    return father_geom_->Global(transform_ * local +
                                offset_.replicate(1, local.cols()));
  }
  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override {
    auto jac = father_geom_->Jacobian(transform_ * local +
                                      offset_.replicate(1, local.cols()));
    for (int i = 0; i < local.cols(); ++i) {
      jac.block(0, 2 * i, DimGlobal(), 2) *= transform_;
    }
    return jac;
  }
  [[nodiscard]] Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override {
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
  [[nodiscard]] Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override {
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
  [[nodiscard]] std::unique_ptr<Geometry> SubGeometry(dim_t codim,
                                                      dim_t i) const override {
    if (codim == 2) {
      return std::make_unique<geometry::Point>(
          Global(base::RefEl::kTria().NodeCoords().col(i)));
    }
    if (codim == 1) {
      Eigen::Matrix2d nodes_in_father;
      nodes_in_father.col(0) = base::RefEl::kTria().NodeCoords().col(
          base::RefEl::kTria().SubSubEntity2SubEntity(1, i, 1, 0));
      nodes_in_father.col(1) = base::RefEl::kTria().NodeCoords().col(
          base::RefEl::kTria().SubSubEntity2SubEntity(1, i, 1, 1));
      return std::make_unique<SegmentIn2DFather>(
          std::make_unique<TriaInFather>(*this), nodes_in_father);
    }
    LF_ASSERT_MSG(codim == 0, "");
    LF_ASSERT_MSG(i == 0, "");
    return std::make_unique<TriaInFather>(*this);
  }
  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const geometry::RefinementPattern& ref_pat,
      lf::base::dim_t codim) const override {
    std::vector<std::unique_ptr<Geometry>> result;
    result.reserve(ref_pat.NumChildren(codim));
    double lattice_const = ref_pat.LatticeConst();
    if (codim == 2) {
      for (auto& n : ref_pat.ChildPolygons(2)) {
        LF_ASSERT_MSG(n.cols() == 1, "");
        result.push_back(std::make_unique<geometry::Point>(
            Global(n.cast<double>() / lattice_const)));
      }
    } else if (codim == 1) {
      for (auto& e : ref_pat.ChildPolygons(1)) {
        LF_ASSERT_MSG(e.cols() == 2, "");
        result.push_back(std::make_unique<SegmentIn2DFather>(
            father_geom_->SubGeometry(0, 0),
            transform_ * e.cast<double>() / lattice_const +
                offset_.replicate<1, 2>()));
      }
    } else if (codim == 0) {
      for (auto& e : ref_pat.ChildPolygons(0)) {
        LF_ASSERT_MSG(e.cols() == 3, "Only triangles are supported so far.");
        result.push_back(std::make_unique<TriaInFather>(
            father_geom_->SubGeometry(0, 0),
            transform_ * e.cast<double>() / lattice_const +
                offset_.replicate<1, 3>()));
      }
    }
    return result;
  }

 private:
  std::unique_ptr<geometry::Geometry> father_geom_;
  Eigen::Matrix<double, 2, 2> transform_;
  Eigen::Vector2d offset_;
};

BrepTriaTransfinite::BrepTriaTransfinite(
    std::array<std::pair<std::shared_ptr<const interface::BrepGeometry>,
                         Eigen::RowVector2d>,
               3>
        curves)
    : curves_(curves) {
  LF_ASSERT_MSG(
      curves_[0].first->DimGlobal() == curves_[1].first->DimGlobal() &&
          curves_[0].first->DimGlobal() == curves_[2].first->DimGlobal(),
      "Inconsistent global dimension.");
  LF_ASSERT_MSG(curves_[0].first->DimLocal() == 1, "curves[0] is not a curve.");
  LF_ASSERT_MSG(curves_[1].first->DimLocal() == 1, "curves[0] is not a curve.");
  LF_ASSERT_MSG(curves_[2].first->DimLocal() == 1, "curves[0] is not a curve.");

  /*LF_ASSERT_MSG(
      curves_[0].first != curves_[1].first ||
          curves_[1].first != curves_[2].first,
      "All three boundaing curves of this triangle are the same curve. This is "
      "not allowed because it yields a degenerate triangle.");*/

  // Make sure the curves agree at the endpoints:
  LF_ASSERT_MSG(
      curves_[0]
          .first->Global(curves_[0].second.col(1))
          .isApprox(curves_[1].first->Global(curves_[1].second.col(0))),
      "curves don't match in node 1: \n"
          << curves_[0].first->Global(curves_[0].second.col(1)).transpose()
          << " <-> "
          << curves_[1].first->Global(curves_[1].second.col(0)).transpose());
  LF_ASSERT_MSG(
      curves_[1]
          .first->Global(curves_[1].second.col(1))
          .isApprox(curves_[2].first->Global(curves_[2].second.col(0))),
      "curves don't match in node 2:\n"
          << curves_[1].first->Global(curves_[1].second.col(1)).transpose()
          << " <-> "
          << curves_[2].first->Global(curves_[2].second.col(0)).transpose());
  LF_ASSERT_MSG(
      curves_[2]
          .first->Global(curves_[2].second.col(1))
          .isApprox(curves_[0].first->Global(curves_[0].second.col(0))),
      "curves don't match in node 0:\n"
          << curves_[2].first->Global(curves_[2].second.col(1)).transpose()
          << " <-> "
          << curves_[0].first->Global(curves_[0].second.col(0)).transpose());
  node0_ = curves_[0].first->Global(curves_[0].second.col(0));
}

Eigen::MatrixXd BrepTriaTransfinite::Global(
    const Eigen::MatrixXd& local) const {
  auto g0 = [&](auto& x) {
    return curves_[0]
        .first
        ->Global((curves_[0].second(0) +
                  (curves_[0].second(1) - curves_[0].second(0)) * x)
                     .matrix())
        .array()
        .eval();
  };
  auto g1 = [&](auto& x) {
    return curves_[1]
        .first
        ->Global((curves_[1].second(0) +
                  (curves_[1].second(1) - curves_[1].second(0)) * x)
                     .matrix())
        .array()
        .eval();
  };
  auto g2 = [&](auto& x) {
    return curves_[2]
        .first
        ->Global((curves_[2].second(0) +
                  (curves_[2].second(1) - curves_[2].second(0)) * x)
                     .matrix())
        .array()
        .eval();
  };
  auto x = local.row(0).array().eval();
  auto y = local.row(1).array().eval();
  auto xr = x.replicate(DimGlobal(), 1).eval();
  auto yr = y.replicate(DimGlobal(), 1).eval();

  return (g0(x) * (1 - yr) - xr * g0(1 - y) +
          (node0_ * (x + y - 1).matrix()).array() + g2(1 - y) * (1 - xr) -
          yr * g2(x) + xr * g1(y) + yr * g1(1 - x))
      .matrix();
}

Eigen::MatrixXd BrepTriaTransfinite::Jacobian(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 2, "expected 2d coordinate pairs");
  Eigen::MatrixXd result(DimGlobal(), 2 * local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    auto g0 = [&](double x) {
      return curves_[0].first->Global(
          curves_[0].second.col(0) +
          (curves_[0].second.col(1) - curves_[0].second.col(0)) * x);
    };
    auto g1 = [&](double x) {
      return curves_[1].first->Global(
          curves_[1].second.col(0) +
          (curves_[1].second.col(1) - curves_[1].second.col(0)) * x);
    };
    auto g2 = [&](double x) {
      return curves_[2].first->Global(
          curves_[2].second.col(0) +
          (curves_[2].second.col(1) - curves_[2].second.col(0)) * x);
    };
    auto Dg0 = [&](double x) -> Eigen::VectorXd {
      return curves_[0].first->Jacobian(
                 curves_[0].second.col(0) +
                 (curves_[0].second.col(1) - curves_[0].second.col(0)) * x) *
             (curves_[0].second(1) - curves_[0].second(0));
    };
    auto Dg1 = [&](double x) -> Eigen::VectorXd {
      return curves_[1].first->Jacobian(
                 curves_[1].second.col(0) +
                 (curves_[1].second.col(1) - curves_[1].second.col(0)) * x) *
             (curves_[1].second.col(1) - curves_[1].second.col(0));
    };
    auto Dg2 = [&](double x) -> Eigen::VectorXd {
      return curves_[2].first->Jacobian(
                 curves_[2].second.col(0) +
                 (curves_[2].second.col(1) - curves_[2].second.col(0)) * x) *
             (curves_[2].second.col(1) - curves_[2].second.col(0));
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
    return std::make_unique<BrepCurve>(curves_[i].first, curves_[i].second);
  }

  // codim==2
  LF_ASSERT_MSG(i >= 0 && i < 3, "subidx out of range.");
  if (i == 0) {
    return std::make_unique<geometry::Point>(node0_);
  } else {
    return std::make_unique<geometry::Point>(
        curves_[i].first->Global(curves_[i].second.col(0)));
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

      result.push_back(std::make_unique<TriaInFather>(
          std::make_unique<BrepTriaTransfinite>(*this),
          c.cast<double>() / lattice_const));
    }
  } else if (codim == 1) {
    for (auto& c : ref_pat.ChildPolygons(1)) {
      LF_ASSERT_MSG(c.cols() == 2, "Expected 2 columns...");
      result.push_back(std::make_unique<SegmentIn2DFather>(
          std::make_unique<BrepTriaTransfinite>(*this),
          c.cast<double>() / lattice_const));
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
BrepTriaTransfinitePerronnet::BrepTriaTransfinitePerronnet(
    std::array<std::pair<std::shared_ptr<const interface::BrepGeometry>,
                         Eigen::RowVector2d>,
               3>
        curves)
    : curves_(curves), nodes_(curves_[0].first->DimGlobal(), 3) {
  LF_ASSERT_MSG(
      curves_[0].first->DimGlobal() == curves_[1].first->DimGlobal() &&
          curves_[0].first->DimGlobal() == curves_[2].first->DimGlobal(),
      "Inconsistent global dimension.");
  LF_ASSERT_MSG(curves_[0].first->DimLocal() == 1, "curves[0] is not a curve.");
  LF_ASSERT_MSG(curves_[1].first->DimLocal() == 1, "curves[0] is not a curve.");
  LF_ASSERT_MSG(curves_[2].first->DimLocal() == 1, "curves[0] is not a curve.");

  /*LF_ASSERT_MSG(
      curves_[0].first != curves_[1].first ||
          curves_[1].first != curves_[2].first,
      "All three boundaing curves of this triangle are the same curve. This is "
      "not allowed because it yields a degenerate triangle.");*/

  // Make sure the curves agree at the endpoints:
  LF_ASSERT_MSG(
      (curves_[0].first->Global(curves_[0].second.col(1)) -
       curves_[1].first->Global(curves_[1].second.col(0)))
              .norm() < 1e-6,
      "curves don't match in node 1: \n"
          << curves_[0].first->Global(curves_[0].second.col(1)).transpose()
          << " <-> "
          << curves_[1].first->Global(curves_[1].second.col(0)).transpose());
  LF_ASSERT_MSG(
      (curves_[1].first->Global(curves_[1].second.col(1)) -
       curves_[2].first->Global(curves_[2].second.col(0)))
              .norm() < 1e-6,
      "curves don't match in node 2:\n"
          << curves_[1].first->Global(curves_[1].second.col(1)).transpose()
          << " <-> "
          << curves_[2].first->Global(curves_[2].second.col(0)).transpose());
  LF_ASSERT_MSG(
      (curves_[2].first->Global(curves_[2].second.col(1)) -
       curves_[0].first->Global(curves_[0].second.col(0)))
              .norm() < 1e-6,
      "curves don't match in node 0:\n"
          << curves_[2].first->Global(curves_[2].second.col(1)).transpose()
          << " <-> "
          << curves_[0].first->Global(curves_[0].second.col(0)).transpose());

  nodes_.col(0) = curves_[0].first->Global(curves_[0].second.col(0));
  nodes_.col(1) = curves_[1].first->Global(curves_[1].second.col(0));
  nodes_.col(2) = curves_[2].first->Global(curves_[2].second.col(0));
}

Eigen::MatrixXd BrepTriaTransfinitePerronnet::Global(
    const Eigen::MatrixXd& local) const {
  auto g0 = [&](auto& x) {
    return curves_[0]
        .first
        ->Global((curves_[0].second(0) +
                  (curves_[0].second(1) - curves_[0].second(0)) * x)
                     .matrix())
        .array()
        .eval();
  };
  auto g1 = [&](auto& x) {
    return curves_[1]
        .first
        ->Global((curves_[1].second(0) +
                  (curves_[1].second(1) - curves_[1].second(0)) * x)
                     .matrix())
        .array()
        .eval();
  };
  auto g2 = [&](auto& x) {
    return curves_[2]
        .first
        ->Global((curves_[2].second(0) +
                  (curves_[2].second(1) - curves_[2].second(0)) * x)
                     .matrix())
        .array()
        .eval();
  };
  auto x = local.row(0).array().eval();
  auto y = local.row(1).array().eval();
  auto xr = x.replicate(DimGlobal(), 1).eval();
  auto yr = y.replicate(DimGlobal(), 1).eval();

  return ((1 - xr - yr) * (g0(x) + g2(1 - y) -
                           nodes_.col(0).replicate(1, local.cols()).array()) +
          xr * (g1(y) + g0(x + y) -
                nodes_.col(1).replicate(1, local.cols()).array()) +
          yr * (g2(1 - x - y) + g1(1 - x) -
                nodes_.col(2).replicate(1, local.cols()).array()))
      .matrix();
}

Eigen::MatrixXd BrepTriaTransfinitePerronnet::Jacobian(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 2, "expected 2d coordinate pairs");
  Eigen::MatrixXd result(DimGlobal(), 2 * local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    auto g0 = [&](double x) {
      return curves_[0].first->Global(
          curves_[0].second.col(0) +
          (curves_[0].second.col(1) - curves_[0].second.col(0)) * x);
    };
    auto g1 = [&](double x) {
      return curves_[1].first->Global(
          curves_[1].second.col(0) +
          (curves_[1].second.col(1) - curves_[1].second.col(0)) * x);
    };
    auto g2 = [&](double x) {
      return curves_[2].first->Global(
          curves_[2].second.col(0) +
          (curves_[2].second.col(1) - curves_[2].second.col(0)) * x);
    };
    auto Dg0 = [&](double x) -> Eigen::VectorXd {
      return curves_[0].first->Jacobian(
                 curves_[0].second.col(0) +
                 (curves_[0].second.col(1) - curves_[0].second.col(0)) * x) *
             (curves_[0].second(1) - curves_[0].second(0));
    };
    auto Dg1 = [&](double x) -> Eigen::VectorXd {
      return curves_[1].first->Jacobian(
                 curves_[1].second.col(0) +
                 (curves_[1].second.col(1) - curves_[1].second.col(0)) * x) *
             (curves_[1].second.col(1) - curves_[1].second.col(0));
    };
    auto Dg2 = [&](double x) -> Eigen::VectorXd {
      return curves_[2].first->Jacobian(
                 curves_[2].second.col(0) +
                 (curves_[2].second.col(1) - curves_[2].second.col(0)) * x) *
             (curves_[2].second.col(1) - curves_[2].second.col(0));
    };
    auto x = local(0, i);
    auto y = local(1, i);

    result.col(2 * i) = nodes_.col(0) - nodes_.col(1) - g0(x) + g0(x + y) +
                        g1(y) - g2(1 - y) - (x + y - 1) * Dg0(x) +
                        x * Dg0(x + y) - y * (Dg1(1 - x) + Dg2(1 - x - y));
    result.col(2 * i + 1) = nodes_.col(0) - nodes_.col(2) - g0(x) + g1(1 - x) -
                            g2(1 - y) + g2(1 - x - y) +
                            x * (Dg0(x + y) + Dg1(y)) +
                            (x + y - 1) * Dg2(1 - y) - y * Dg2(1 - x - y);
  }
  return result;
}

Eigen::MatrixXd BrepTriaTransfinitePerronnet::JacobianInverseGramian(
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

Eigen::VectorXd BrepTriaTransfinitePerronnet::IntegrationElement(
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

std::unique_ptr<geometry::Geometry> BrepTriaTransfinitePerronnet::SubGeometry(
    dim_t codim, dim_t i) const {
  LF_ASSERT_MSG(codim >= 0 && codim <= 2, "codim out of range.");
  if (codim == 0) {
    LF_ASSERT_MSG(i == 0, "subidx out of range.");
    return std::make_unique<BrepTriaTransfinitePerronnet>(*this);
  } else if (codim == 1) {
    LF_ASSERT_MSG(i >= 0 && i < 3, "subidx out of range.");
    return std::make_unique<BrepCurve>(curves_[i].first, curves_[i].second);
  }

  // codim==2
  LF_ASSERT_MSG(i >= 0 && i < 3, "subidx out of range.");
  return std::make_unique<geometry::Point>(nodes_.col(i));
}

std::vector<std::unique_ptr<geometry::Geometry>>
BrepTriaTransfinitePerronnet::ChildGeometry(
    const geometry::RefinementPattern& ref_pat, lf::base::dim_t codim) const {
  LF_ASSERT_MSG(codim >= 0 && codim < 3, "codim out of bounds.");
  std::vector<std::unique_ptr<geometry::Geometry>> result;
  result.reserve(ref_pat.NumChildren(codim));
  double lattice_const = static_cast<double>(ref_pat.LatticeConst());
  if (codim == 0) {
    for (auto& c : ref_pat.ChildPolygons(0)) {
      LF_ASSERT_MSG(
          c.cols() == 3,
          "So far a triangle can only be split into other triangles...");

      result.push_back(std::make_unique<TriaInFather>(
          std::make_unique<BrepTriaTransfinitePerronnet>(*this),
          c.cast<double>() / lattice_const));
    }
  } else if (codim == 1) {
    for (auto& c : ref_pat.ChildPolygons(1)) {
      LF_ASSERT_MSG(c.cols() == 2, "Expected 2 columns...");
      result.push_back(std::make_unique<SegmentIn2DFather>(
          std::make_unique<BrepTriaTransfinitePerronnet>(*this),
          c.cast<double>() / lattice_const));
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
}  // namespace lf::brep::geom

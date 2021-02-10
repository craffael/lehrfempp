/**
 * @file
 * @brief Definition of BrepLagrSegment
 * @author Raffael Casagrande
 * @date   2021-01-28 05:05:01
 * @copyright MIT License
 */

#ifndef __e4bffaf701e04343919102118e3df626
#define __e4bffaf701e04343919102118e3df626

#include <lf/brep/interface/interface.h>
#include <lf/geometry/geometry.h>

namespace lf::brep::geom {

template <class GEOMETRY>
class BrepLagrSegment : public geometry::Geometry {
 public:
  BrepLagrSegment(GEOMETRY segment,
                  std::shared_ptr<const interface::BrepCurve> curve,
                  double curve_coord_p0, double curve_coord_p1);
  BrepLagrSegment(const BrepLagrSegment& other) = default;

  [[nodiscard]] dim_t DimLocal() const override { return geom_.DimLocal(); }
  [[nodiscard]] dim_t DimGlobal() const override { return geom_.DimGlobal(); }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kSegment();
  }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override {
    return geom_.Global(local);
  }
  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override {
    return geom_.Jacobian(local);
  }
  [[nodiscard]] Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override {
    return geom_.JacobianInverseGramian(local);
  }
  [[nodiscard]] Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override {
    return geom_.IntegrationElement(local);
  }
  [[nodiscard]] std::unique_ptr<Geometry> SubGeometry(dim_t codim,
                                                      dim_t i) const override {
    LF_ASSERT_MSG(codim >= 0 && codim < 2, "Codim out of range.");
    if (codim == 0) {
      LF_ASSERT_MSG(i == 0, "index i out of range.");
      return std::make_unique<BrepLagrSegment>(*this);
    }
    return geom_.SubGeometry(codim, i);
  }

  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const geometry::RefinementPattern& ref_pat,
      lf::base::dim_t codim) const override;

  static const auto LagrangeNodes() { return GEOMETRY::LagrangeNodes(); }

 private:
  GEOMETRY geom_;
  std::shared_ptr<const brep::interface::BrepCurve> curve_;
  double curve_coord_p0_;
  double curve_coord_p1_;
};

template <class GEOMETRY>
BrepLagrSegment<GEOMETRY>::BrepLagrSegment(
    GEOMETRY segment, std::shared_ptr<const interface::BrepCurve> curve,
    double curve_coord_p0,

    double curve_coord_p1)
    : geom_(std::move(segment)),
      curve_(curve),
      curve_coord_p0_(curve_coord_p0),
      curve_coord_p1_(curve_coord_p1) {
  LF_ASSERT_MSG(geom_.RefEl() == base::RefEl::kSegment(),
                "BrepLagrSegment can only wrap a geometry of type segment.");
  LF_ASSERT_MSG((geom_.LagrangeNodes().array() >= 0.).all() &&
                    (geom_.LagrangeNodes().array() <= 1).all(),
                "lagr_nodes are out of bounds.");
  LF_ASSERT_MSG(std::abs(geom_.LagrangeNodes()(0)) < 1e-6,
                "lagr_nodes(0) node must coincide with left segment end.");
  LF_ASSERT_MSG(std::abs(geom_.LagrangeNodes()(1) - 1) < 1e-6,
                "lagr_nodes(1) must coincide with right segment end.");
  LF_ASSERT_MSG(
      curve->GlobalSingle(curve_coord_p0_)
          .isApprox(geom_.Global((Eigen::Matrix<double, 1, 1>::Constant(0.)))),
      "curve(curve_coord_p0) != geom(0)");
  LF_ASSERT_MSG(
      curve->GlobalSingle(curve_coord_p1_)
          .isApprox(geom_.Global((Eigen::Matrix<double, 1, 1>::Constant(1.)))),
      "curve(curve_coord_p1) != geom(1)");
}

template <class GEOMETRY>
std::vector<std::unique_ptr<geometry::Geometry>>
BrepLagrSegment<GEOMETRY>::ChildGeometry(
    const geometry::RefinementPattern& ref_pat, lf::base::dim_t codim) const {
  auto children = ref_pat.ChildPolygons(codim);
  std::vector<std::unique_ptr<geometry::Geometry>> result;
  result.reserve(children.size());
  auto lattice_const = ref_pat.LatticeConst();
  if (codim == 0) {
    for (auto& child : children) {
      LF_ASSERT_MSG(child.cols() == 2, "Expected exactly two points.")
      LF_ASSERT_MSG(child.rows() == 1, "Expected exactly one row.");

      auto curve_local =
          ((geom_.LagrangeNodes().array() * (child(1) - child(0)) + child(0)) /
               static_cast<double>(lattice_const) *
               (curve_coord_p1_ - curve_coord_p0_) +
           curve_coord_p0_)
              .eval();
      result.push_back(std::make_unique<BrepLagrSegment<GEOMETRY>>(
          GEOMETRY(curve_->GlobalMulti(curve_local.matrix())), curve_,
          curve_local(0), curve_local(1)));
    }
  } else {
    LF_ASSERT_MSG(codim == 1, "codim out of bounds.");
    for (auto& child : children) {
      LF_ASSERT_MSG(child.cols() == 1, "Expected exactly one point.");
      LF_ASSERT_MSG(child.rows() == 1, "Expected exactly one row.");
      double curve_local =
          (curve_coord_p1_ - curve_coord_p0_) / lattice_const * child(0) +
          curve_coord_p0_;
      result.push_back(
          std::make_unique<geometry::Point>(curve_->GlobalSingle(curve_local)));
    }
  }
  return result;
}

}  // namespace lf::brep::geom

#endif  // __e4bffaf701e04343919102118e3df626

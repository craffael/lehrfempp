#include "point.h"

namespace lf::geometry {
Eigen::MatrixXd Point::Global(const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 0, "local.rows() != 0");
  return coord_.replicate(1, local.cols());
}

// NOLINTNEXTLINE(misc-unused-parameters)
Eigen::MatrixXd Point::Jacobian(const Eigen::MatrixXd& local) const {
  return Eigen::MatrixXd::Zero(DimGlobal(), 0);
}

Eigen::MatrixXd Point::JacobianInverseGramian(
    const ::Eigen::MatrixXd& local) const {  // NOLINT(misc-unused-parameters)
  return Eigen::MatrixXd::Zero(DimGlobal(), 0);
}

// NOLINTNEXTLINE(misc-unused-parameters)
Eigen::VectorXd Point::IntegrationElement(const Eigen::MatrixXd& local) const {
  return Eigen::Matrix<double,1,1>::Constant(1.0);
}

std::unique_ptr<Geometry> Point::SubGeometry(dim_t codim, dim_t i) const {
  if (codim == 0 && i == 0) {
    return std::make_unique<Point>(coord_);
  }
  LF_VERIFY_MSG(false, "codim or i out of bounds.");
}

std::vector<std::unique_ptr<Geometry>>
Point::ChildGeometry(const RefinementPattern& ref_pattern,
		     lf::base::dim_t codim) const {
  LF_VERIFY_MSG(codim == 0, "Only codim = 0 allowed for a point");
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
  LF_VERIFY_MSG(ref_pattern.RefEl() == lf::base::RefEl::kPoint(),
                "ref_patern.RefEl() = " << ref_pattern.RefEl().ToString());
  LF_VERIFY_MSG(ref_pattern.noChildren(0) == 1,
                "ref_pattern.noChildren() = " << ref_pattern.noChildren(0));
  // The only way to "refine" a point is to copy it
  child_geo_uptrs.push_back(std::make_unique<Point>(coord_));
  return child_geo_uptrs;
}

}  // namespace lf::geometry

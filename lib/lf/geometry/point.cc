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
  LF_VERIFY_MSG(false, "integrationElement() undefined for points.");
}

std::unique_ptr<Geometry> Point::SubGeometry(dim_t codim, dim_t i) const {
  if (codim == 0 && i == 0) {
    return std::make_unique<Point>(coord_);
  }
  LF_VERIFY_MSG(false, "codim or i out of bounds.");
}

  std::vector<std::unique_ptr<Geometry>>
  Point::ChildGeometry(const RefinementPattern &ref_pattern) const
  {
    std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
    LF_VERIFY_MSG(ref_pattern.refpat() == (int)RefPat::rp_copy,
		  "Illegal refinement pattern " << (int)ref_pattern.refpat());
    child_geo_uptrs.push_back(std::make_unique<Point>(coord_));
    return std::move(child_geo_uptrs);
  }
  
  
}  // namespace lf::geometry
